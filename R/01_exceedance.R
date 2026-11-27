# Conteo de excedencias por celda y métricas horarias para un archivo.
# Excedencia estricta: LOAD > umbral. NA/relleno = ausente (no cuenta como cero).
# Equivalente a python/exceedance.py (Accumulator, process_file).

fmt_time <- function(t) format(t, "%Y-%m-%d %H:%M")

# Acumulador de conteos por celda y métricas horarias en una rejilla (lat, lon).
# thrs: umbrales redondeados a float32; thr_names: umbrales originales (para nombrar columnas)
new_accumulator <- function(lat, lon, thrs, thr_names) {
  e <- new.env()
  e$lat <- lat; e$lon <- lon; e$nlon <- length(lon); e$nlat <- length(lat)
  e$area <- cell_area_km2(lat, lon); e$area_vec <- as.numeric(e$area)
  e$counts <- lapply(thrs, function(t) matrix(0L, e$nlon, e$nlat))
  e$valid <- matrix(0L, e$nlon, e$nlat)
  e$rows <- list()
  e$update <- function(z, hour_index) {                 # z: [lon, lat, n]
    n <- dim(z)[3]
    zm <- matrix(z, nrow = e$nlon * e$nlat, ncol = n)    # celdas x horas
    e$valid <- e$valid + matrix(as.integer(rowSums(!is.na(zm))), e$nlon, e$nlat)
    row <- data.table(hour_index = hour_index)
    for (i in seq_along(thrs)) {
      exc <- zm > thrs[i]; exc[is.na(exc)] <- FALSE
      e$counts[[i]] <- e$counts[[i]] + matrix(rowSums(exc), e$nlon, e$nlat)
      s <- thr_suffix(thr_names[i], thr_names[1])
      row[[paste0("ncells", s)]] <- as.integer(colSums(exc))
      row[[paste0("area_km2", s)]] <- as.numeric(crossprod(exc, e$area_vec))
    }
    zi <- zm; zi[is.na(zi)] <- -Inf
    idx <- apply(zi, 2, which.max)                      # celda del máximo por hora
    mx <- zi[cbind(idx, seq_len(n))]; mx[!is.finite(mx)] <- NA
    row$max_load_kg_m2 <- mx
    row$max_lat <- e$lat[(idx - 1) %/% e$nlon + 1]; row$max_lon <- e$lon[(idx - 1) %% e$nlon + 1]
    e$rows[[length(e$rows) + 1]] <- row
  }
  e$result <- function(scenario, year, rec_ok) {
    df <- rbindlist(e$rows)
    t0 <- as.POSIXct(sprintf("%d-01-01", year), tz = "UTC")
    df[, `:=`(scenario = scenario, year = year)]
    df[, time := t0 + hour_index * 3600]
    df[, `:=`(month = as.integer(format(time, "%m")), hour_of_day = as.integer(format(time, "%H")))]
    df[, time := format(time, "%Y-%m-%d %H:%M:%S")]
    df <- df[rec_ok]                                    # registros sin ningún dato no son pronósticos
    setcolorder(df, c("scenario", "year", "hour_index", "time", "max_load_kg_m2", "max_lat", "max_lon"))  # orden de Python
    list(counts = e$counts[[1]], valid = e$valid, hourly = df, lat = e$lat, lon = e$lon, area = e$area)
  }
  e
}

process_file <- function(path, scenario, year, cfg) {
  # LOAD está en float32: se compara contra el umbral redondeado a float32, de modo que un
  # valor almacenado igual a 0.1f no cuente como excedencia (igual que en Python).
  thrs <- to_float32(cfg$thresholds)
  f <- open_load(path, cfg$variable_candidates); on.exit(nc_close(f$nc))
  ver <- detect_version(f$attrs, f$var, cfg$version_override[[as.character(year)]])
  hr <- resolve_hours(f$horas, f$nrec, year)
  warns <- hr$warns
  if (ver$label == "desconocida")
    warns <- c(warns, "no se pudo determinar la versión de Fall3D (usa version_override en config.yml)")
  dlat <- cfg$domain_lat; dlon <- cfg$domain_lon
  if (abs(min(f$lat) - dlat[1]) > 0.01 || abs(max(f$lat) - dlat[2]) > 0.01 ||
      abs(min(f$lon) - dlon[1]) > 0.01 || abs(max(f$lon) - dlon[2]) > 0.01)
    warns <- c(warns, "el dominio del archivo difiere del configurado (domain_lat/domain_lon)")

  # Rejilla de trabajo (y nativa, si se interpola: los combinados de una sola versión se
  # calculan sin interpolar)
  rg <- NULL
  if (cfg$grid_policy == "regrid") {
    tg <- target_axes(cfg)
    if (!is_same_grid(f$lat, f$lon, tg$lat, tg$lon)) rg <- make_regridder(f$lat, f$lon, tg$lat, tg$lon)
  }
  acc_nat <- new_accumulator(f$lat, f$lon, thrs, cfg$thresholds)
  acc_work <- if (!is.null(rg)) new_accumulator(tg$lat, tg$lon, thrs, cfg$thresholds) else acc_nat
  grid_native <- sprintf("%dx%d", f$nlat, f$nlon)
  grid_used <- sprintf("%dx%d", acc_work$nlat, acc_work$nlon)
  nh <- f$nrec
  rec_ok <- logical()
  inv <- list(n_nan = 0, n_negative = 0, vmin = Inf, vmax = -Inf)

  for (h0 in seq(1, nh, by = cfg$chunk_hours)) {
    n <- min(cfg$chunk_hours, nh - h0 + 1)
    z <- read_load_chunk(f, h0, n, cfg$fill_limit)      # NaN, Inf, relleno o enmascarado -> NA
    ok <- !is.na(z)
    inv$n_nan <- inv$n_nan + sum(!ok)
    inv$n_negative <- inv$n_negative + sum(z < 0, na.rm = TRUE)
    if (any(ok)) { inv$vmin <- min(inv$vmin, min(z, na.rm = TRUE)); inv$vmax <- max(inv$vmax, max(z, na.rm = TRUE)) }
    rec_ok <- c(rec_ok, colSums(matrix(ok, nrow = f$nlon * f$nlat, ncol = n)) > 0)  # registro con datos (nativa)
    hidx <- hr$idx[h0:(h0 + n - 1)]
    acc_nat$update(z, hidx)
    if (!is.null(rg)) acc_work$update(regrid_block(z, rg), hidx)   # interpolar LOAD antes del umbral
  }
  work <- acc_work$result(scenario, year, rec_ok)
  native <- if (!is.null(rg)) acc_nat$result(scenario, year, rec_ok) else work

  # Estado del año a partir de las horas con datos (no se fija a mano)
  t0 <- as.POSIXct(sprintf("%d-01-01", year), tz = "UTC")
  exp <- expected_hours(year)
  n_valid <- sum(rec_ok)
  if (n_valid == 0) {
    status <- "sin datos"; last_hour <- -1L; first_hour <- -1L
  } else {
    last_hour <- max(hr$idx[rec_ok]); first_hour <- min(hr$idx[rec_ok])
    status <- if (last_hour < exp - 1) "parcial" else if (n_valid < exp) "con huecos" else "completo"
  }
  n_missing <- if (n_valid > 0) last_hour + 1L - n_valid else 0L
  last_time <- if (n_valid > 0) fmt_time(t0 + last_hour * 3600) else ""
  first_time <- if (n_valid > 0) fmt_time(t0 + first_hour * 3600) else ""
  parts <- character()
  if (status == "parcial") parts <- c(parts, sprintf("parcial (hasta %s)", last_time))
  if (status == "sin datos") parts <- c(parts, "sin datos")
  if (n_missing > 0) parts <- c(parts, sprintf("%d h sin datos", n_missing))
  if (paste0(scenario, "-", year) %in% cfg$replaced_files && nzchar(cfg$replaced_note))
    parts <- c(parts, cfg$replaced_note)     # registros reemplazados (ver docs/correccion_c0XM.md)
  inv <- c(list(scenario = scenario, year = year, file = basename(path), variable = f$var,
                fall3d_version = ver$label, fall3d_version_full = ver$full, version_source = ver$source,
                grid_native = grid_native, grid = grid_used, regridded = !is.null(rg),
                n_hours = nh, n_valid_hours = as.integer(n_valid), expected_full_year = exp,
                first_time = first_time, last_time = last_time, time_source = hr$source,
                n_horas_mismatch = hr$n_mismatch, n_missing_hours = as.integer(n_missing), status = status,
                complete = status == "completo", partial_flag = status == "parcial",
                nota = paste(parts, collapse = "; ")), inv)
  for (w in warns) message(sprintf("[AVISO] %s: %s", basename(path), w))
  list(work = work, native = native, inv = as.data.table(inv))
}
