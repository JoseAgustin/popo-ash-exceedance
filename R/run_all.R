#!/usr/bin/env Rscript
# Punto de entrada.
# Autor: Dr. José Agustín García Reynoso (ICAyCC, UNAM) <agustin@atmosfera.unam.mx>
# Licencia: MIT (ver LICENSE).
# Uso (desde la raíz del repositorio):
#   Rscript R/run_all.R [config/config.yml] [output_dir]
# data_dir se toma de la variable de entorno POPO_DATA_DIR si está definida; si no, de config.yml.
# Equivalente a python/exceedance.py: cada archivo se inspecciona por separado (variable LOAD o
# tephra_grn_load, rejilla, versión de Fall3D, tiempo) y grid_policy decide cómo tratar 100x100 vs 101x101.
args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[1] else "config/config.yml"

for (f in c("00_utils.R", "01_exceedance.R", "02_metrics.R", "03_maps_plots.R", "04_export.R"))
  source(file.path("R", f))

cfg <- read_config(cfg_path)
if (nzchar(Sys.getenv("POPO_DATA_DIR"))) cfg$data_dir <- Sys.getenv("POPO_DATA_DIR")
if (!dir.exists(cfg$data_dir))
  stop("data_dir no existe: '", cfg$data_dir, "'. Edita config/config.yml o define la variable de entorno POPO_DATA_DIR.")
outdir <- if (length(args) >= 2) args[2] else cfg$output_dir
for (d in c("tables", "netcdf")) dir.create(file.path(outdir, d), recursive = TRUE, showWarnings = FALSE)

# Dos espacios de datos: "work" (rejilla de trabajo) y "native" (rejilla de cada archivo, sin interpolar)
SP <- c("work", "native")
P <- setNames(lapply(SP, function(sp) list(counts = list(), valid = list())), SP)
G <- setNames(lapply(SP, function(sp) list()), SP)
hourly <- setNames(lapply(SP, function(sp) list()), SP)
inv <- list()
t0 <- Sys.time()
for (s in cfg$scenarios) for (y in cfg$years) {
  path <- file.path(cfg$data_dir, gsub("{year}", y, gsub("{scenario}", s, cfg$file_pattern, fixed = TRUE), fixed = TRUE))
  if (!file.exists(path)) { message(sprintf("[AVISO] no existe %s; se omite", basename(path))); next }
  r <- process_file(path, s, y, cfg)
  k <- paste(s, y, sep = "|")
  for (sp in SP) {
    P[[sp]]$counts[[k]] <- r[[sp]]$counts; P[[sp]]$valid[[k]] <- r[[sp]]$valid
    hourly[[sp]][[k]] <- r[[sp]]$hourly
    if (is.null(G[[sp]][[as.character(y)]])) G[[sp]][[as.character(y)]] <- list(lat = r[[sp]]$lat, lon = r[[sp]]$lon, area = r[[sp]]$area)
  }
  inv[[k]] <- r$inv
  i <- r$inv
  message(sprintf("[%5.0fs] %s: Fall3D %s, %s, %s->%s, %d registros, %s, horas con exc = %d, max = %.3f kg/m2",
                  as.numeric(difftime(Sys.time(), t0, units = "secs")), basename(path), i$fall3d_version, i$variable,
                  i$grid_native, i$grid, i$n_hours, i$status, sum(r$work$hourly$ncells > 0), i$vmax))
}
if (!length(inv)) stop("No se procesó ningún archivo.")

inv_dt <- rbindlist(inv)
thrs <- cfg$thresholds
# Información por año (se verifica que los escenarios coincidan en versión y rejilla)
yinfo <- list()
for (y in sort(unique(inv_dt$year))) {
  g <- inv_dt[year == y]
  for (col in c("fall3d_version", "grid", "grid_native"))
    if (uniqueN(g[[col]]) > 1) message(sprintf("[AVISO] %d: los escenarios difieren en %s: %s", y, col, paste(g[[col]], collapse = ", ")))
  worst <- if (any(g$status != "completo")) g[status != "completo"][1] else g[1]
  yinfo[[as.character(y)]] <- list(version = g$fall3d_version[1], grid = g$grid[1], grid_native = g$grid_native[1],
                                   regridded = g$regridded[1], status = worst$status, last_time = worst$last_time,
                                   n_missing = worst$n_missing_hours)
}
combos <- build_combos(yinfo)
for (cmb in combos) if (!cmb$ok)
  message(sprintf("[AVISO] periodo %s (%s) omitido: sus años usan rejillas distintas (%s); usa grid_policy: regrid",
                  cmb$label, cmb$kind, cmb$grid))

# Información (versión, rejilla, nota) por escenario y periodo
info <- inv_dt[, .(scenario, period = as.character(year), fall3d_version, grid, grid_native, nota)]
ystat <- lapply(cfg$scenarios, function(s) {
  z <- inv_dt[scenario == s]
  setNames(lapply(seq_len(nrow(z)), function(i) list(status = z$status[i], last_time = z$last_time[i],
                                                       n_missing = z$n_missing_hours[i])), as.character(z$year))
}); names(ystat) <- cfg$scenarios
extra <- list()
for (cmb in combos) if (cmb$ok) for (s in cfg$scenarios) {
  if (length(ystat[[s]]) == 0) next
  extra[[length(extra) + 1]] <- data.table(scenario = s, period = cmb$label, fall3d_version = cmb$version, grid = cmb$grid,
    grid_native = cmb$grid_native, nota = combo_note(cmb, as.integer(names(ystat[[s]])), ystat[[s]]))
}
info <- rbindlist(c(list(info), extra))

H <- lapply(hourly, rbindlist); Hp <- expand_periods(H, combos)
by_year <- attach_info(summary_table(Hp, c("scenario", "period"), thrs), info, c("scenario", "period"))
monthly <- attach_info(summary_table(Hp, c("scenario", "period", "month"), thrs), info, c("scenario", "period", "month"))
Hp[, season := ifelse(month %in% cfg$wet_months, "lluvias(may-oct)", "secas(nov-abr)")]
seasonal <- attach_info(summary_table(Hp, c("scenario", "period", "season"), thrs), info, c("scenario", "period", "season"))
diurnal <- attach_info(summary_table(Hp, c("scenario", "period", "hour_of_day"), thrs), info, c("scenario", "period", "hour_of_day"))
pt <- prob_tables(P, G, cfg, combos)
max_prob <- attach_info(pt$max_prob, info, c("scenario", "period"))
area_levels <- attach_info(pt$area_levels, info, c("scenario", "period", "prob_level"))

pdefs <- rbindlist(c(
  lapply(names(yinfo), function(y) data.table(period = y, kind = "año", years = y, fall3d_version = yinfo[[y]]$version,
         grid = yinfo[[y]]$grid, grid_native = yinfo[[y]]$grid_native, status = yinfo[[y]]$status)),
  lapply(combos, function(cmb) data.table(period = cmb$label, kind = cmb$kind, years = paste(cmb$years, collapse = " "),
         fall3d_version = cmb$version, grid = cmb$grid, grid_native = cmb$grid_native,
         status = if (cmb$ok) "combinado" else "omitido (rejillas distintas)"))))

Hw <- H$work
td <- file.path(outdir, "tables")
fwrite(inv_dt, file.path(td, "file_inventory.csv"))
fwrite(pdefs, file.path(td, "period_definitions.csv"))
fwrite(by_year, file.path(td, "summary_scenario_period.csv"))
fwrite(monthly, file.path(td, "summary_monthly.csv"))
fwrite(seasonal, file.path(td, "summary_seasonal.csv"))
fwrite(diurnal, file.path(td, "summary_diurnal.csv"))
fwrite(max_prob, file.path(td, "max_probability.csv"))
fwrite(area_levels, file.path(td, "area_by_probability_level.csv"))
fwrite(Hw[ncells > 0, .(scenario, year, hour_index, time, ncells, area_km2, max_load_kg_m2, max_lat, max_lon)],
       file.path(td, "exceedance_hours.csv"))
fwrite(Hw, file.path(td, "hourly_metrics.csv.gz"))

for (p in netcdf_plan(yinfo, combos)) {
  desc <- if (p$space == "native") "native grid (not interpolated)" else
    if (cfg$grid_policy == "regrid") "common grid (regrid)" else "native grid"
  g <- G[[p$space]][[as.character(p$years[1])]]
  write_prob_netcdf(file.path(outdir, "netcdf", p$name), cfg, g$lat, g$lon, P[[p$space]], p$years, yinfo, p$sets, desc)
}

valid_combos <- Filter(function(c) c$ok, combos)
primary <- if (length(valid_combos)) valid_combos[[1]]$label else names(yinfo)[1]
make_plots(outdir, cfg, G, P, yinfo, combos,
           list(by_year = by_year, monthly = monthly[period == primary], diurnal = diurnal[period == primary], hourly = Hw))
message("Listo. Salidas en ", normalizePath(outdir))
