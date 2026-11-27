# Utilidades: configuración, lectura robusta de NetCDF (Fall3D 7.2 y v8), versión, tiempo y rejilla.
# Autor: Dr. José Agustín García Reynoso (ICAyCC, UNAM) <agustin@atmosfera.unam.mx>. Licencia: MIT.
suppressPackageStartupMessages({
  library(ncdf4)
  library(yaml)
  library(data.table)
})

EARTH_RADIUS_KM <- 6371.0088

DEFAULTS <- list(
  file_pattern = "{scenario}-{year}.nc",
  variable_candidates = c("LOAD", "tephra_grn_load"),
  version_override = list(),
  grid_policy = "regrid",
  regrid_n = 100,
  domain_lon = c(-100.62, -96.62),
  domain_lat = c(17.02, 21.02),
  extra_thresholds = list()
)

read_config <- function(path = "config/config.yml") {
  cfg <- yaml::read_yaml(path)
  for (k in names(DEFAULTS)) if (is.null(cfg[[k]])) cfg[k] <- list(DEFAULTS[[k]])
  cfg$extra_thresholds <- as.numeric(unlist(cfg$extra_thresholds))
  cfg$fill_limit <- as.numeric(cfg$fill_limit)  # yaml puede leer 1.0e30 como texto
  cfg$thresholds <- c(as.numeric(cfg$threshold), cfg$extra_thresholds)
  cfg$variable_candidates <- as.character(unlist(cfg$variable_candidates))
  cfg$domain_lon <- as.numeric(unlist(cfg$domain_lon)); cfg$domain_lat <- as.numeric(unlist(cfg$domain_lat))
  cfg$years <- as.integer(unlist(cfg$years))
  if (!cfg$grid_policy %in% c("native", "regrid"))
    stop("grid_policy debe ser 'native' o 'regrid', no '", cfg$grid_policy, "'")
  ov <- list()
  for (item in unlist(cfg$version_override)) {
    p <- strsplit(as.character(item), ":", fixed = TRUE)[[1]]
    ov[[p[1]]] <- trimws(p[2])
  }
  cfg$version_override <- ov
  cfg
}

# Redondea números a precisión float32 (devuelve double); vectorizado
to_float32 <- function(x) {
  x <- as.numeric(x)
  readBin(writeBin(x, raw(), size = 4), "numeric", size = 4, n = length(x))
}

# Sufijo de columnas para umbrales adicionales (el principal no lleva sufijo)
thr_suffix <- function(t, main) if (t == main) "" else paste0("_gt", format(t, trim = TRUE))

# Área (km2) de cada celda; devuelve matriz [lon, lat] (orden de ncdf4)
cell_area_km2 <- function(lat, lon) {
  dlat <- abs(mean(diff(lat))) * pi / 180
  dlon <- abs(mean(diff(lon))) * pi / 180
  a_lat <- EARTH_RADIUS_KM^2 * dlat * dlon * cos(lat * pi / 180)
  matrix(rep(a_lat, each = length(lon)), nrow = length(lon), ncol = length(lat))
}

# ---------------------------------------------------------------------------
# Versión de Fall3D y tiempo de cada registro
# ---------------------------------------------------------------------------
version_label <- function(full) {
  p <- strsplit(full, ".", fixed = TRUE)[[1]]
  if (as.integer(p[1]) >= 8) p[1] else paste(p[1], p[2], sep = ".")
}

# Devuelve list(label, full, source). Etiqueta: "7.2", "8" o "desconocida"
detect_version <- function(attrs, varname, override = NULL) {
  if (!is.null(override)) return(list(label = if (startsWith(override, "8")) "8" else override,
                                      full = override, source = "config:version_override"))
  grab <- function(key) if (!is.null(attrs[[key]])) as.character(attrs[[key]]) else NULL
  title <- grab("TITLE")
  if (!is.null(title)) {
    m <- regmatches(title, regexec("fall3d\\s*v?([0-9]+\\.[0-9]+(\\.[0-9]+)?)", title, ignore.case = TRUE))[[1]]
    if (length(m)) return(list(label = version_label(m[2]), full = m[2], source = "TITLE"))
  }
  for (key in c("Code_version", "code_version", "source")) {
    v <- grab(key)
    if (!is.null(v)) {
      m <- regmatches(v, regexec("([0-9]+\\.[0-9]+(\\.[0-9]+)?)", v))[[1]]
      if (length(m)) return(list(label = version_label(m[2]), full = m[2], source = key))
    }
  }
  if (varname == "tephra_grn_load") return(list(label = "8", full = "8", source = "variable"))
  list(label = "desconocida", full = "desconocida", source = "ninguno")
}

expected_hours <- function(year) if ((year %% 4 == 0 && year %% 100 != 0) || year %% 400 == 0) 8784L else 8760L

# Hora del año (desde 0) de cada registro. list(idx, source, n_mismatch, warns)
resolve_hours <- function(horas, nrec, year) {
  idx <- 0:(nrec - 1)
  warns <- character()
  if (is.null(horas)) return(list(idx = idx, source = "index", n_mismatch = 0L, warns = warns))
  h <- as.numeric(horas)
  if (length(h) <= 1 || (length(h) == nrec && all(h == h[1])))
    return(list(idx = idx, source = "index", n_mismatch = 0L, warns = warns))   # escalar o constante
  if (length(h) != nrec)
    return(list(idx = idx, source = "index", n_mismatch = 0L,
                warns = sprintf("'horas' tiene %d valores y hay %d registros; se usa el índice", length(h), nrec)))
  hi <- round(h)
  if (any(abs(h - hi) > 1e-6) || any(diff(hi) <= 0) || hi[1] < 0 || hi[length(hi)] >= expected_hours(year))
    return(list(idx = idx, source = "index", n_mismatch = 0L,
                warns = "'horas' no es entera, creciente y dentro del año; se usa el índice"))
  n_mis <- sum(hi != idx)
  if (n_mis > 0) warns <- sprintf("'horas' difiere del índice en %d registros (huecos en la serie); se usa 'horas'", n_mis)
  list(idx = hi, source = "horas", n_mismatch = as.integer(n_mis), warns = warns)
}

# ---------------------------------------------------------------------------
# Apertura y lectura por bloques
# ---------------------------------------------------------------------------
# Abre el archivo y detecta variable, dimensiones, rejilla, atributos y 'horas'.
open_load <- function(path, var_candidates) {
  nc <- nc_open(path)
  vn <- var_candidates[var_candidates %in% names(nc$var)]
  if (!length(vn)) {
    nc_close(nc)
    stop(basename(path), ": no hay ninguna de las variables ", paste(var_candidates, collapse = ", "),
         "; hay ", paste(names(nc$var), collapse = ", "))
  }
  vn <- vn[1]
  dims <- nc$var[[vn]]$dim                       # orden de R: inverso al de NetCDF
  dn <- vapply(dims, function(d) d$name, "")
  dl <- tolower(dn)
  ilat <- which(startsWith(dl, "lat"))[1]; ilon <- which(startsWith(dl, "lon"))[1]
  irec <- length(dims)                           # la primera dimensión de NetCDF es el registro
  if (is.na(ilat) || is.na(ilon) || irec %in% c(ilat, ilon))
    stop(basename(path), ": ", vn, " debe tener dimensiones lat/lon y el registro como primera dimensión")
  coord <- function(i) as.numeric(ncvar_get(nc, dn[i]))
  list(nc = nc, var = vn, dims = dims, ilat = ilat, ilon = ilon, irec = irec,
       nrec = dims[[irec]]$len, nlat = dims[[ilat]]$len, nlon = dims[[ilon]]$len,
       lat = coord(ilat), lon = coord(ilon),
       attrs = ncatt_get(nc, 0),
       horas = if ("horas" %in% names(nc$var)) as.numeric(ncvar_get(nc, "horas")) else NULL)
}

# Lee un bloque de registros: array [lon, lat, n] en double. Dimensiones extra (time) -> último valor.
read_load_chunk <- function(f, h0, n, fill_limit) {
  nd <- length(f$dims)
  start <- integer(nd); count <- integer(nd)
  for (i in seq_len(nd)) {
    len <- f$dims[[i]]$len
    if (i == f$irec) { start[i] <- h0; count[i] <- n }
    else if (i %in% c(f$ilat, f$ilon)) { start[i] <- 1L; count[i] <- len }
    else { start[i] <- len; count[i] <- 1L }
  }
  z <- ncvar_get(f$nc, f$var, start = start, count = count, collapse_degen = FALSE)
  others <- setdiff(seq_len(nd), c(f$ilon, f$ilat, f$irec))
  z <- array(aperm(z, c(f$ilon, f$ilat, f$irec, others)), dim = c(f$nlon, f$nlat, n))
  z[!is.finite(z) | z >= fill_limit] <- NA_real_
  z
}

# ---------------------------------------------------------------------------
# Rejilla común (interpolación bilineal)
# ---------------------------------------------------------------------------
target_axes <- function(cfg) {
  n <- cfg$regrid_n
  list(lat = seq(cfg$domain_lat[1], cfg$domain_lat[2], length.out = n),
       lon = seq(cfg$domain_lon[1], cfg$domain_lon[2], length.out = n))
}

is_same_grid <- function(lat, lon, tlat, tlon, tol = 1e-4) {
  length(lat) == length(tlat) && length(lon) == length(tlon) &&
    max(abs(lat - tlat)) < tol && max(abs(lon - tlon)) < tol
}

# Índices (base 1) del nodo inferior, peso del superior y marca "fuera del dominio"
axis_weights <- function(src, tgt, tol = 1e-4) {
  if (any(diff(src) <= 0)) stop("las coordenadas de la rejilla deben ser estrictamente crecientes")
  i0 <- pmin(pmax(findInterval(tgt, src), 1L), length(src) - 1L)
  w <- pmin(pmax((tgt - src[i0]) / (src[i0 + 1L] - src[i0]), 0), 1)
  list(i0 = i0, w = w, outside = tgt < src[1] - tol | tgt > src[length(src)] + tol)
}

make_regridder <- function(src_lat, src_lon, tgt_lat, tgt_lon) {
  y <- axis_weights(src_lat, tgt_lat); x <- axis_weights(src_lon, tgt_lon)
  list(iy = y$i0, wy = y$w, ix = x$i0, wx = x$w,
       outside = outer(x$outside, y$outside, `|`), nlon_t = length(tgt_lon), nlat_t = length(tgt_lat))
}

# z: array [lon, lat, n] (double). Interpola en doble precisión (primero lat, luego lon) y redondea a float32.
regrid_block <- function(z, rg) {
  wy <- rep(rg$wy, each = dim(z)[1])                       # [lon, lat_t], se recicla sobre n
  b <- (1 - wy) * z[, rg$iy, , drop = FALSE] + wy * z[, rg$iy + 1L, , drop = FALSE]
  wx <- rep(rg$wx, times = rg$nlat_t)                      # [lon_t, lat_t]
  cc <- (1 - wx) * b[rg$ix, , , drop = FALSE] + wx * b[rg$ix + 1L, , , drop = FALSE]
  cc[rep(rg$outside, dim(cc)[3])] <- NA_real_
  array(to_float32(cc), dim = dim(cc))
}

# Mapa base (estados CONABIO, costa Natural Earth, recortado) desde data/basemap/*.csv
load_basemap <- function(dir = "data/basemap") {
  rd <- function(k) { f <- file.path(dir, paste0(k, ".csv")); if (file.exists(f)) fread(f) else NULL }
  list(states = rd("states"), coast = rd("coast"), state_labels = rd("state_labels"), cities = rd("cities"))
}
