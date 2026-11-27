# Mapas y gráficos (ggplot2). Equivalente a python/exceedance.py (make_plots).
# Mapa base: estados CONABIO/INEGI y costa Natural Earth (data/basemap/*.csv).
suppressPackageStartupMessages(library(ggplot2))

# Probabilidad por celda en formato largo; w/h = tamaño de celda (cada panel puede tener su rejilla)
prob_df <- function(counts, valid, lat, lon, label) {
  d <- expand.grid(lon = lon, lat = lat)
  d$prob <- as.numeric(ifelse(valid > 0, counts / pmax(valid, 1), NA))
  d$w <- mean(diff(lon)); d$h <- mean(diff(lat))
  d$panel <- label
  d
}

# Capas de mapa base + ciudades + cráter + límites del dominio
basemap_layers <- function(bm, cfg, lon, lat, labels = TRUE) {
  L <- list()
  if (!is.null(bm$states)) L <- c(L, list(geom_path(data = bm$states, aes(lon, lat, group = id),
                                                     colour = "grey25", linewidth = 0.3, inherit.aes = FALSE)))
  if (!is.null(bm$coast)) L <- c(L, list(geom_path(data = bm$coast, aes(lon, lat, group = id),
                                                    colour = "black", linewidth = 0.6, inherit.aes = FALSE)))
  if (labels && !is.null(bm$state_labels))
    L <- c(L, list(geom_text(data = bm$state_labels, aes(lon, lat, label = name), size = 2.2,
                             fontface = "italic", colour = "grey35", inherit.aes = FALSE)))
  if (labels && !is.null(bm$cities)) {
    ci <- copy(bm$cities)
    ci[, `:=`(tx = lon + ifelse(ha == "left", 0.05, -0.05), ty = lat + ifelse(dy < 0, -0.04, 0.04),
              hj = ifelse(ha == "left", 0, 1), vj = ifelse(dy < 0, 1, 0))]
    L <- c(L, list(geom_point(data = ci, aes(lon, lat), shape = 21, fill = "white", size = 1.3,
                              inherit.aes = FALSE),
                   geom_text(data = ci, aes(tx, ty, label = name, hjust = hj, vjust = vj), size = 2.5,
                             inherit.aes = FALSE)))
  }
  c(L, list(geom_point(data = data.frame(lon = cfg$crater_lon, lat = cfg$crater_lat), aes(lon, lat),
                       shape = 24, fill = "blue", colour = "black", size = 2.6, inherit.aes = FALSE),
            coord_fixed(ratio = 1 / cos(19 * pi / 180), xlim = range(lon), ylim = range(lat), expand = FALSE)))
}

partial_text <- function(ys, yinfo) {
  txt <- character()
  for (y in ys) {
    st <- yinfo[[as.character(y)]]
    if (st$status == "parcial") txt <- c(txt, sprintf("%d parcial (hasta %s)", y, st$last_time))
    else if (st$status == "con huecos") txt <- c(txt, sprintf("%d con huecos (%d h sin datos)", y, st$n_missing))
    else if (st$status != "completo") txt <- c(txt, sprintf("%d %s", y, st$status))
  }
  paste(txt, collapse = "; ")
}

make_plots <- function(outdir, cfg, G, P, yinfo, combos, tabs) {
  fig <- file.path(outdir, "figures"); dir.create(fig, recursive = TRUE, showWarnings = FALSE)
  thr <- cfg$threshold
  years <- sort(as.integer(names(yinfo)))
  heights <- setNames(cfg$column_height_km, cfg$scenarios)
  bm <- load_basemap()
  foot <- "Estados: CONABIO/INEGI; costa: Natural Earth"
  rng <- function(sp, y) list(lon = G[[sp]][[as.character(y)]]$lon, lat = G[[sp]][[as.character(y)]]$lat)
  first_version <- yinfo[[as.character(years[1])]]$version

  valid_combos <- Filter(function(c) c$ok, combos)
  maps <- if (length(valid_combos)) valid_combos else
    list(list(kind = "todos", years = years[1], label = as.character(years[1]), version = first_version,
              space = "work", tag = "all", grid = yinfo[[as.character(years[1])]]$grid))
  has_homog <- any(vapply(maps, function(c) c$kind == "homogeneo", TRUE))

  combo_df <- function(cmb) {
    g <- rng(cmb$space, cmb$years[1]); Ps <- P[[cmb$space]]
    rbindlist(lapply(cfg$scenarios, function(s) {
      kk <- paste(s, cmb$years, sep = "|"); kk <- kk[kk %in% names(Ps$counts)]
      if (!length(kk)) return(NULL)
      prob_df(Reduce(`+`, Ps$counts[kk]), Reduce(`+`, Ps$valid[kk]), g$lat, g$lon,
              sprintf("%s (columna %d km), %s", s, heights[[s]], cmb$label))
    }))
  }

  for (cmb in maps) {
    fname <- if (cmb$kind == "homogeneo") {
      if (cmb$version == first_version) "maps_prob_combined.png" else sprintf("maps_prob_combined_%s.png", cmb$tag)
    } else if (has_homog) "maps_prob_combined_all.png" else "maps_prob_combined.png"
    d <- combo_df(cmb); if (is.null(d) || !nrow(d)) next
    d$prob[d$prob == 0] <- NA
    g <- rng(cmb$space, cmb$years[1])
    sub <- sprintf("Fall3D %s, rejilla %s%s%s", cmb$version, cmb$grid, if (cmb$kind == "homogeneo") " nativa" else "",
                   if (grepl("+", cmb$version, fixed = TRUE)) " (mezcla de versiones: puede reflejar el cambio del modelo)" else "")
    pt <- partial_text(cmb$years, yinfo)
    p <- ggplot(d, aes(lon, lat, fill = prob, width = w, height = h)) + geom_tile() +
      scale_fill_distiller(palette = "YlOrRd", direction = 1, trans = "sqrt", na.value = "white",
                           name = sprintf("P(LOAD > %g kg/m2)", thr)) +
      basemap_layers(bm, cfg, g$lon, g$lat) + facet_wrap(~panel) + theme_bw() +
      labs(title = sprintf("Probabilidad de excedencia (fracción de horas). Triángulo = cráter. %s", foot),
           subtitle = paste0(sub, if (nzchar(pt)) paste0(". ", pt) else ""))
    ggsave(file.path(fig, fname), p, width = 17, height = 6.4, dpi = 150)
  }

  Pw <- P$work
  for (s in cfg$scenarios) {
    ys <- Filter(function(y) paste(s, y, sep = "|") %in% names(Pw$counts), years)
    if (!length(ys)) next
    d <- rbindlist(lapply(ys, function(y) {
      kk <- paste(s, y, sep = "|"); g <- rng("work", y); st <- yinfo[[as.character(y)]]
      extra <- if (st$status == "parcial") sprintf("\nparcial, hasta %s", st$last_time)
               else if (st$status == "con huecos") sprintf("\ncon huecos (%d h sin datos)", st$n_missing)
               else if (st$status != "completo") paste0("\n", st$status) else ""
      lab <- sprintf("%s %d (Fall3D %s%s)%s", s, y, st$version, if (isTRUE(st$regridded)) ", interpolado" else "", extra)
      prob_df(Pw$counts[[kk]], Pw$valid[[kk]], g$lat, g$lon, lab)
    }))
    d$prob[d$prob == 0] <- NA
    g <- rng("work", ys[1])
    p <- ggplot(d, aes(lon, lat, fill = prob, width = w, height = h)) + geom_tile() +
      scale_fill_distiller(palette = "YlOrRd", direction = 1, trans = "sqrt", na.value = "white",
                           name = sprintf("P(LOAD > %g)", thr)) +
      basemap_layers(bm, cfg, g$lon, g$lat, labels = FALSE) + facet_wrap(~panel, ncol = 4) + theme_bw()
    ggsave(file.path(fig, sprintf("maps_prob_by_year_%s.png", s)), p, width = 16, height = 4 * ceiling(length(ys) / 4) + 1, dpi = 130)
  }

  # Diferencias entre escenarios (primer periodo combinado válido)
  c0 <- maps[[1]]; g0 <- rng(c0$space, c0$years[1]); P0 <- P[[c0$space]]
  pcomb <- lapply(cfg$scenarios, function(s) {
    kk <- paste(s, c0$years, sep = "|"); kk <- kk[kk %in% names(P0$counts)]
    if (!length(kk)) return(NULL)
    prob_df(Reduce(`+`, P0$counts[kk]), Reduce(`+`, P0$valid[kk]), g0$lat, g0$lon, s)$prob
  }); names(pcomb) <- cfg$scenarios
  pairs <- Filter(function(p) !is.null(pcomb[[p[1]]]) && !is.null(pcomb[[p[2]]]),
                  list(c("c05M", "c03M"), c("c10M", "c05M"), c("c10M", "c03M")))
  if (length(pairs)) {
    grid <- expand.grid(lon = g0$lon, lat = g0$lat)
    dd <- rbindlist(lapply(pairs, function(p) cbind(grid, diff = pcomb[[p[1]]] - pcomb[[p[2]]],
                                                     w = mean(diff(g0$lon)), h = mean(diff(g0$lat)),
                                                     panel = paste(p[1], "-", p[2]))))
    dm <- max(abs(dd$diff), na.rm = TRUE)
    p <- ggplot(dd, aes(lon, lat, fill = diff, width = w, height = h)) + geom_tile() +
      scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-dm, dm), name = "Diferencia de\nprobabilidad") +
      basemap_layers(bm, cfg, g0$lon, g0$lat) + facet_wrap(~panel) + theme_bw() +
      labs(title = sprintf("Diferencia entre escenarios, %s (Fall3D %s). %s", c0$label, c0$version, foot))
    ggsave(file.path(fig, "maps_prob_difference.png"), p, width = 17, height = 6.4, dpi = 150)
  }

  by_year <- tabs$by_year[period %in% as.character(years)]
  by_year[, estado := vapply(period, function(y) yinfo[[y]]$status, "")]
  by_year[, etiqueta := vapply(period, function(y) sprintf("%s%s\nFall3D %s", y, if (yinfo[[y]]$status != "completo") "*" else "", yinfo[[y]]$version), "")]
  by_year[, etiqueta := factor(etiqueta, levels = unique(etiqueta[order(period)]))]
  pt <- partial_text(years, yinfo)
  p <- ggplot(by_year, aes(etiqueta, 100 * frac_hours_exc, fill = scenario, alpha = estado == "completo")) +
    geom_col(position = "dodge") + scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.55), guide = "none") +
    theme_bw() + labs(x = "Año (* = parcial o con huecos)", y = sprintf("%% de horas con LOAD > %g en alguna celda", thr),
                      caption = pt)
  ggsave(file.path(fig, "hours_exceedance_by_year.png"), p, width = max(9, 1.2 * length(years) + 3), height = 4.8, dpi = 150)

  H <- tabs$hourly[ncells > 0]
  p <- ggplot(H, aes(scenario, area_km2)) + geom_boxplot() + scale_y_log10() + theme_bw() +
    labs(y = "Área con LOAD > umbral (km2), horas con excedencia",
         title = "Área de influencia por escenario (todos los años, rejilla de trabajo)")
  ggsave(file.path(fig, "area_by_scenario.png"), p, width = 6.5, height = 4.5, dpi = 150)

  ttl <- sprintf("%s (Fall3D %s)", c0$label, c0$version)
  p <- ggplot(tabs$monthly, aes(month, 100 * frac_hours_exc, colour = scenario)) +
    geom_line() + geom_point() + theme_bw() + labs(x = "Mes", y = "% de horas con excedencia", title = ttl)
  ggsave(file.path(fig, "monthly_exceedance.png"), p, width = 7.5, height = 4.2, dpi = 150)
  p <- ggplot(tabs$diurnal, aes(hour_of_day, 100 * frac_hours_exc, colour = scenario)) +
    geom_line() + geom_point() + theme_bw() + labs(x = "Hora del día", y = "% de horas con excedencia", title = ttl)
  ggsave(file.path(fig, "diurnal_exceedance.png"), p, width = 7.5, height = 4.2, dpi = 150)
}
