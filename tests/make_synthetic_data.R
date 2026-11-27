#!/usr/bin/env Rscript
# Crea archivos NetCDF sintéticos pequeños que imitan Fall3D 7.2 y v8, y sus config.yml, FUERA del repo.
# Autor: Dr. José Agustín García Reynoso (ICAyCC, UNAM) <agustin@atmosfera.unam.mx>. Licencia: MIT.
# Uso: Rscript tests/make_synthetic_data.R <carpeta_destino>
#
# dataset1 (rejillas reales 100x100 y 101x101, pocas horas, config por defecto):
#   2021, 2022: Fall3D 7.2, LOAD, 100x100, 'horas' escalar 0, 72 registros
#   2024: v8 (atributo source 8.0.1), variable tephra_grn_load, 101x101, 70 registros con 'horas'
#         con un hueco de 4 h (0:39, 44:73) y un registro todo relleno
#   2025: v8 (Code_version 8.2.1), variable LOAD, 101x101, 72 registros, 'horas' = índice
# dataset2 (rejillas 10x10 y 11x11, años completos, regrid_n = 10):
#   2023: 7.2, 8760 registros (completo). 2024: v8, 8784 registros (bisiesto, completo).
#   2025: v8, 5000 registros (parcial)
# nc4: un archivo NetCDF-4 (solo para probar R; Python clásico no lo lee sin netCDF4)
suppressPackageStartupMessages(library(ncdf4))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Uso: Rscript tests/make_synthetic_data.R <carpeta_destino>")
root <- args[1]

# Campo continuo común: lo muestrean todas las rejillas
field <- function(lat, lon, h) {
  amp <- 1.6 * abs(sin(0.37 * h + 1))
  x0 <- -98.6 + 0.5 * sin(0.11 * h); y0 <- 19.0 + 0.4 * cos(0.13 * h)
  outer(lon, lat, function(x, y) amp * exp(-(((x - x0) / 0.5)^2 + ((y - y0) / 0.4)^2))) - 0.002
}

write_nc <- function(path, n, year, nrec, varname = "LOAD", style = c("7.2", "v8_source", "v8_code"),
                     horas = c("scalar0", "index", "custom"), hvals = NULL, nan_rec = integer(), force_v4 = FALSE) {
  style <- match.arg(style); horas <- match.arg(horas)
  lat <- seq(17.02, 21.02, length.out = n); lon <- seq(-100.62, -96.62, length.out = n)
  d_lon <- ncdim_def("lon", "degrees_east", lon); d_lat <- ncdim_def("lat", "degrees_north", lat)
  d_time <- ncdim_def("time", sprintf("hours since %d-01-01 00:00:0.0", year), 0)
  d_hora <- ncdim_def("hora", "", 1:nrec, unlim = TRUE, create_dimvar = FALSE)
  v_load <- ncvar_def(varname, "kg/m2", list(d_lon, d_lat, d_time, d_hora), prec = "float")
  v_horas <- if (horas == "scalar0") ncvar_def("horas", "", list(), prec = "integer") else
    ncvar_def("horas", "", list(d_hora), prec = "integer")
  nc <- nc_create(path, list(v_load, v_horas), force_v4 = force_v4)
  hv <- switch(horas, scalar0 = 0L, index = 0:(nrec - 1), custom = as.integer(hvals))
  hh <- if (horas == "scalar0") 0:(nrec - 1) else hv           # hora real de cada registro (campo)
  arr <- array(0, c(n, n, 1, nrec))
  for (i in seq_len(nrec)) arr[, , 1, i] <- field(lat, lon, hh[i])
  for (i in nan_rec) arr[, , 1, i] <- NA
  ncvar_put(nc, varname, arr, start = c(1, 1, 1, 1), count = c(n, n, 1, nrec))   # dim ilimitada: count explícito
  if (horas == "scalar0") ncvar_put(nc, "horas", hv) else ncvar_put(nc, "horas", hv, start = 1, count = nrec)
  if (style == "7.2") ncatt_put(nc, 0, "TITLE", "Fall3d 7.2 results")
  if (style == "v8_source") ncatt_put(nc, 0, "source", "FALL3D model version 8.0.1")
  if (style == "v8_code") ncatt_put(nc, 0, "Code_version", "8.2.1")
  nc_close(nc)
}

cfg_text <- function(data_dir, out_dir, years, n_regrid, grid_policy = "regrid", chunk = 50) {
  c(sprintf("data_dir: %s", data_dir), sprintf("output_dir: %s", out_dir),
    "scenarios: [c03M]", "column_height_km: [3]", sprintf("years: [%s]", paste(years, collapse = ", ")),
    "file_pattern: \"{scenario}-{year}.nc\"", "variable_candidates: [LOAD, tephra_grn_load]", "version_override: []",
    sprintf("grid_policy: %s", grid_policy), sprintf("regrid_n: %d", n_regrid),
    "domain_lon: [-100.62, -96.62]", "domain_lat: [17.02, 21.02]",
    "threshold: 0.1", "extra_thresholds: [1.0]", "fill_limit: 1.0e+30", sprintf("chunk_hours: %d", chunk),
    "crater_lat: 19.023", "crater_lon: -98.622", "prob_levels: [0.001, 0.01, 0.05]", "wet_months: [5, 6, 7, 8, 9, 10]")
}

d1 <- file.path(root, "dataset1"); dir.create(d1, recursive = TRUE, showWarnings = FALSE)
write_nc(file.path(d1, "c03M-2021.nc"), 100, 2021, 72, "LOAD", "7.2", "scalar0")
write_nc(file.path(d1, "c03M-2022.nc"), 100, 2022, 72, "LOAD", "7.2", "scalar0")
write_nc(file.path(d1, "c03M-2024.nc"), 101, 2024, 70, "tephra_grn_load", "v8_source", "custom",
         hvals = c(0:39, 44:73), nan_rec = 11)
write_nc(file.path(d1, "c03M-2025.nc"), 101, 2025, 72, "LOAD", "v8_code", "index")
file.create(file.path(d1, "c03M-2024.nc.prev_20240424"))      # copia que el patrón no debe recoger
writeLines(cfg_text(d1, file.path(root, "out_d1_regrid"), c(2021, 2022, 2023, 2024, 2025), 100),
           file.path(root, "dataset1_regrid.yml"))
writeLines(cfg_text(d1, file.path(root, "out_d1_native"), c(2021, 2022, 2023, 2024, 2025), 100, "native"),
           file.path(root, "dataset1_native.yml"))

d2 <- file.path(root, "dataset2"); dir.create(d2, recursive = TRUE, showWarnings = FALSE)
write_nc(file.path(d2, "c03M-2023.nc"), 10, 2023, 8760, "LOAD", "7.2", "scalar0")
write_nc(file.path(d2, "c03M-2024.nc"), 11, 2024, 8784, "tephra_grn_load", "v8_source", "index")
write_nc(file.path(d2, "c03M-2025.nc"), 11, 2025, 5000, "LOAD", "v8_code", "index")
writeLines(cfg_text(d2, file.path(root, "out_d2"), c(2023, 2024, 2025), 10, chunk = 1000),
           file.path(root, "dataset2.yml"))

d3 <- file.path(root, "nc4"); dir.create(d3, recursive = TRUE, showWarnings = FALSE)
write_nc(file.path(d3, "c03M-2025.nc"), 101, 2025, 72, "LOAD", "v8_code", "index", force_v4 = TRUE)
writeLines(cfg_text(d3, file.path(root, "out_nc4"), 2025, 100), file.path(root, "nc4.yml"))
cat("Datos sintéticos creados en", root, "\n")
