# Resúmenes por escenario/periodo (año y combinados), mes, estación y hora del día.
# Equivalente a python/exceedance.py (build_combos, expand_periods, attach_info, ...).

summarize_group <- function(g, thrs) {
  main <- thrs[1]
  out <- list(n_hours = nrow(g), max_load_kg_m2 = max(g$max_load_kg_m2, na.rm = TRUE))
  for (t in thrs) {
    s <- thr_suffix(t, main)
    nc <- g[[paste0("ncells", s)]]; ar <- g[[paste0("area_km2", s)]]
    e <- nc > 0; ae <- ar[e]
    q <- function(p) if (any(e)) as.numeric(quantile(ae, p, type = 7)) else 0
    out[[paste0("hours_exc", s)]] <- sum(e)
    out[[paste0("frac_hours_exc", s)]] <- mean(e)
    out[[paste0("mean_cells_exc_hours", s)]] <- if (any(e)) mean(nc[e]) else 0
    out[[paste0("area_km2_mean_all_hours", s)]] <- mean(ar)
    out[[paste0("area_km2_mean_exc_hours", s)]] <- if (any(e)) mean(ae) else 0
    out[[paste0("area_km2_p50_exc_hours", s)]] <- q(0.5)
    out[[paste0("area_km2_p90_exc_hours", s)]] <- q(0.9)
    out[[paste0("area_km2_p99_exc_hours", s)]] <- q(0.99)
    out[[paste0("area_km2_max", s)]] <- max(ar)
  }
  out
}

summary_table <- function(H, keys, thrs) {
  H[, summarize_group(.SD, thrs), by = keys, .SDcols = names(H)]
}

# Periodos combinados. Hay dos "espacios" de datos: "work" (rejilla de trabajo) y "native" (rejilla de
# cada archivo, sin interpolar; igual a work en años no interpolados).
#   * combinado homogéneo (una sola versión de Fall3D, p. ej. 2018-2023 (TITLE 7.2, real 7.1) o 2024-2025 con v8):
#     espacio native (requiere que sus años compartan rejilla nativa);
#   * combinado de todos los años (mezcla de versiones): espacio work (requiere rejilla común).
# Con una sola versión solo hay el combinado "todos". yinfo: lista nombrada por año (character).
version_tag <- function(v) paste0("v", gsub(".", "", v, fixed = TRUE))

build_combos <- function(yinfo) {
  years <- sort(as.integer(names(yinfo)))
  if (length(years) < 2) return(list())
  vers <- vapply(as.character(years), function(y) yinfo[[y]]$version, "")
  mk <- function(kind, ys, space) {
    key <- if (space == "native") "grid_native" else "grid"
    grids <- unique(vapply(as.character(ys), function(y) yinfo[[y]][[key]], ""))
    version <- paste(unique(vers[as.character(ys)]), collapse = "+")
    list(kind = kind, years = ys, label = sprintf("%d-%d", ys[1], ys[length(ys)]), version = version, space = space,
         tag = if (kind == "homogeneo") version_tag(version) else "all",
         grid = paste(grids, collapse = "+"),
         grid_native = paste(unique(vapply(as.character(ys), function(y) yinfo[[y]]$grid_native, "")), collapse = "+"),
         ok = length(grids) == 1)
  }
  if (length(unique(vers)) == 1) return(list(mk("todos", years, "work")))
  combos <- list()
  for (v in unique(vers)) {
    ys <- years[vers == v]
    if (length(ys) > 1) combos[[length(combos) + 1]] <- mk("homogeneo", ys, "native")
  }
  c(combos, list(mk("todos", years, "work")))
}

combo_note <- function(cmb, scen_years, ystat) {
  parts <- character()
  if (cmb$kind == "homogeneo") parts <- sprintf("solo Fall3D %s, rejilla nativa %s", cmb$version, cmb$grid)
  else if (grepl("+", cmb$version, fixed = TRUE))
    parts <- sprintf("mezcla Fall3D %s: puede reflejar el cambio de versión del modelo", cmb$version)
  for (y in cmb$years) {
    st <- ystat[[as.character(y)]]
    if (y %in% scen_years && st$status != "completo")
      parts <- c(parts, sprintf("incluye año %s %d%s", st$status, y,
                                if (st$status == "parcial") sprintf(" (hasta %s)", st$last_time)
                                else if (st$status == "con huecos") sprintf(" (%d h sin datos)", st$n_missing) else ""))
  }
  paste(parts, collapse = "; ")
}

# H: lista espacio -> data.table horario
expand_periods <- function(H, combos) {
  parts <- list(copy(H$work)[, period := as.character(year)])
  for (cmb in combos) if (cmb$ok) parts[[length(parts) + 1]] <- copy(H[[cmb$space]][year %in% cmb$years])[, period := cmb$label]
  rbindlist(parts)
}

# Añade versión, rejilla y nota (por escenario y periodo) tras las columnas clave
attach_info <- function(df, info, keys) {
  out <- info[df, on = c("scenario", "period")]
  other <- setdiff(names(df), keys)
  out[, c(keys, "fall3d_version", "grid", "grid_native", other, "nota"), with = FALSE]
}

# Probabilidad máxima y área por nivel de probabilidad (cada año y cada combinado válido)
prob_tables <- function(P, G, cfg, combos) {
  mx <- list(); lv <- list()
  for (s in cfg$scenarios) {
    Pw <- P$work
    ks <- Filter(function(k) startsWith(k, paste0(s, "|")), names(Pw$counts))
    if (!length(ks)) next
    yrs <- sort(as.integer(sub(".*\\|", "", ks)))
    pers <- list()
    for (y in yrs) { k <- paste(s, y, sep = "|"); pers[[as.character(y)]] <- list(c = Pw$counts[[k]], v = Pw$valid[[k]], sp = "work", y0 = y) }
    for (cmb in combos) {
      Ps <- P[[cmb$space]]
      ys <- Filter(function(y) paste(s, y, sep = "|") %in% names(Ps$counts), cmb$years)
      if (cmb$ok && length(ys)) {
        kk <- paste(s, ys, sep = "|")
        pers[[cmb$label]] <- list(c = Reduce(`+`, Ps$counts[kk]), v = Reduce(`+`, Ps$valid[kk]), sp = cmb$space, y0 = ys[1])
      }
    }
    for (p in names(pers)) {
      g <- G[[pers[[p]]$sp]][[as.character(pers[[p]]$y0)]]
      pr <- ifelse(pers[[p]]$v > 0, pers[[p]]$c / pmax(pers[[p]]$v, 1), NA)
      k <- which.max(pr)
      mx[[length(mx) + 1]] <- data.table(scenario = s, period = p, max_prob = max(pr, na.rm = TRUE),
        max_prob_lat = g$lat[(k - 1) %/% length(g$lon) + 1], max_prob_lon = g$lon[(k - 1) %% length(g$lon) + 1],
        cells_ever_exceeded = sum(pers[[p]]$c > 0), area_km2_ever_exceeded = sum(g$area[pers[[p]]$c > 0]))
      for (l in cfg$prob_levels) {
        m <- !is.na(pr) & pr >= l
        lv[[length(lv) + 1]] <- data.table(scenario = s, period = p, prob_level = l,
          n_cells = sum(m), area_km2 = sum(g$area[m]))
      }
    }
  }
  list(max_prob = rbindlist(mx), area_levels = rbindlist(lv))
}
