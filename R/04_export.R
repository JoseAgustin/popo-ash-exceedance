# Exportación a NetCDF. Equivalente a python/exceedance.py (write_netcdf, netcdf_plan).

# NetCDF de un grupo de años con la misma rejilla. sets: lista de list(tag, years, what)
write_prob_netcdf <- function(path, cfg, lat, lon, Pspace, years, yinfo, sets, space_desc) {
  scs <- cfg$scenarios; fv <- -9999
  nlon <- length(lon); nlat <- length(lat)
  ce <- array(0L, c(nlon, nlat, length(years), length(scs))); ve <- ce
  for (i in seq_along(scs)) for (j in seq_along(years)) {
    k <- paste(scs[i], years[j], sep = "|")
    if (!k %in% names(Pspace$counts)) next
    ce[, , j, i] <- Pspace$counts[[k]]; ve[, , j, i] <- Pspace$valid[[k]]
  }
  prob <- function(c, v) ifelse(v > 0, c / pmax(v, 1), fv)

  d_lon <- ncdim_def("lon", "degrees_east", lon); d_lat <- ncdim_def("lat", "degrees_north", lat)
  d_yr <- ncdim_def("year", "", years); d_sc <- ncdim_def("scenario", "", seq_along(scs) - 1L)
  vars <- list(
    ncvar_def("column_height_km", "km", list(d_sc), prec = "integer"),
    ncvar_def("prob_exc_year", "1", list(d_lon, d_lat, d_yr, d_sc), fv, prec = "float",
              longname = "probability of exceedance per year"),
    ncvar_def("count_exc_year", "1", list(d_lon, d_lat, d_yr, d_sc), prec = "integer"),
    ncvar_def("n_valid_year", "1", list(d_lon, d_lat, d_yr, d_sc), prec = "integer"))
  for (st in sets) {
    yl <- paste(st$years, collapse = ",")
    vars <- c(vars, list(
      ncvar_def(paste0("prob_exc_", st$tag), "1", list(d_lon, d_lat, d_sc), fv, prec = "float",
                longname = sprintf("probability of exceedance, %s (%s)", st$what, yl)),
      ncvar_def(paste0("count_exc_", st$tag), "1", list(d_lon, d_lat, d_sc), prec = "integer",
                longname = sprintf("hours with LOAD > threshold, %s (%s)", st$what, yl)),
      ncvar_def(paste0("n_valid_", st$tag), "1", list(d_lon, d_lat, d_sc), prec = "integer",
                longname = sprintf("valid hours, %s (%s)", st$what, yl))))
  }
  nc <- nc_create(path, vars)
  ncvar_put(nc, "column_height_km", as.integer(cfg$column_height_km))
  ncvar_put(nc, "prob_exc_year", prob(ce, ve)); ncvar_put(nc, "count_exc_year", ce); ncvar_put(nc, "n_valid_year", ve)
  for (st in sets) {
    sel <- match(st$years, years)
    ca <- apply(ce[, , sel, , drop = FALSE], c(1, 2, 4), sum); va <- apply(ve[, , sel, , drop = FALSE], c(1, 2, 4), sum)
    ncvar_put(nc, paste0("prob_exc_", st$tag), prob(ca, va))
    ncvar_put(nc, paste0("count_exc_", st$tag), ca); ncvar_put(nc, paste0("n_valid_", st$tag), va)
  }
  ncatt_put(nc, 0, "title", "Probabilidad de excedencia de LOAD (FALL3D, Popocatepetl)")
  ncatt_put(nc, 0, "threshold_kg_m2", cfg$threshold)
  ncatt_put(nc, 0, "scenarios", paste(scs, collapse = ","))
  ncatt_put(nc, 0, "grid", sprintf("%dx%d", nlat, nlon))
  ncatt_put(nc, 0, "grid_policy", cfg$grid_policy)
  ncatt_put(nc, 0, "grid_space", space_desc)
  ncatt_put(nc, 0, "years", paste(years, collapse = ","))
  ncatt_put(nc, 0, "fall3d_version_by_year",
            paste(sprintf("%d:%s", years, vapply(as.character(years), function(y) yinfo[[y]]$version, "")), collapse = ","))
  stt <- vapply(as.character(years), function(y) yinfo[[y]]$status, "")
  pt <- years[stt == "parcial"]
  ncatt_put(nc, 0, "partial_years", if (length(pt)) paste(pt, collapse = ",") else "none")
  gy <- years[stt == "con huecos"]
  ncatt_put(nc, 0, "years_with_gaps", if (length(gy)) paste(sprintf("%d:%dh", gy,
            vapply(as.character(gy), function(y) as.integer(yinfo[[y]]$n_missing), 1L)), collapse = ",") else "none")
  nc_close(nc)
}

# Archivos NetCDF: el de la rejilla de trabajo lleva todos los años (y el combinado de todos); los años
# interpolados tienen además un archivo en su rejilla nativa (con su combinado homogéneo, p. ej. v8).
netcdf_plan <- function(yinfo, combos) {
  years <- sort(as.integer(names(yinfo)))
  yi <- function(y, k) yinfo[[as.character(y)]][[k]]
  plan <- list()
  wg <- split(years, vapply(years, yi, "", k = "grid"))
  wg <- wg[unique(vapply(years, yi, "", k = "grid"))]                 # orden de aparición (como Python)
  for (g in names(wg)) plan[[length(plan) + 1]] <- list(
    name = if (length(wg) == 1) "prob_exceedance.nc" else sprintf("prob_exceedance_%s.nc", g),
    space = "work", grid = g, years = wg[[g]], sets = list())
  ry <- years[vapply(years, yi, TRUE, k = "regridded")]
  if (length(ry)) {
    ng <- split(ry, vapply(ry, yi, "", k = "grid_native"))
    ng <- ng[unique(vapply(ry, yi, "", k = "grid_native"))]
    for (g in names(ng)) plan[[length(plan) + 1]] <- list(name = sprintf("prob_exceedance_%s.nc", g),
      space = "native", grid = g, years = ng[[g]], sets = list())
  }
  for (cmb in combos) {
    if (!cmb$ok) next
    what <- if (cmb$kind == "todos") sprintf("all years (Fall3D %s)", cmb$version) else sprintf("Fall3D %s only, native grid", cmb$version)
    # Un combinado homogéneo sin años interpolados tiene los mismos datos en ambos espacios
    target <- if (cmb$space == "work" || !any(vapply(cmb$years, yi, TRUE, k = "regridded"))) "work" else "native"
    for (i in seq_along(plan)) {
      p <- plan[[i]]
      if (p$space == target && p$grid == cmb$grid && all(cmb$years %in% p$years)) {
        plan[[i]]$sets[[length(p$sets) + 1]] <- list(tag = cmb$tag, years = cmb$years, what = what)
        break
      }
    }
  }
  plan
}
