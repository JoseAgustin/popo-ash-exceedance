#!/usr/bin/env python3
"""Excedencia de depósito de ceniza (LOAD) en salidas FALL3D (7.2 y v8) del Popocatépetl.

Autor: Dr. José Agustín García Reynoso (ICAyCC, UNAM) <agustin@atmosfera.unam.mx>
Licencia: MIT (ver LICENSE).

Lee los NetCDF cXXM-AAAA.nc por bloques de horas, calcula para cada celda la
frecuencia/probabilidad de LOAD > umbral y métricas horarias (celdas excedidas,
área en km2, máximo), y escribe tablas CSV, NetCDF y figuras PNG.

Cada archivo se inspecciona por separado (variable LOAD o tephra_grn_load, dimensiones y
rejilla reales, versión de Fall3D, tiempo de cada registro); ver python/popo_io.py y docs/datos.md.
Rejillas distintas (100x100 en Fall3D 7.2 y 101x101 en v8) se manejan con `grid_policy`
(config.yml): 'regrid' interpola LOAD a una rejilla común ANTES de comparar con el umbral;
'native' deja cada año en su rejilla y solo combina años que comparten rejilla.

Uso:
    python python/exceedance.py --config config/config.yml
    python python/exceedance.py --scenarios c03M --years 2018 --output-dir outputs_test

Lectores NetCDF (en este orden): netCDF4, scipy.io.netcdf_file (solo formato clásico).
"""
import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from popo_io import (NCFile, Regridder, cell_area_km2, detect_version, expected_hours,  # noqa: E402
                     is_same_grid, resolve_hours, target_axes)

DEFAULTS = {
    "file_pattern": "{scenario}-{year}.nc",
    "variable_candidates": ["LOAD", "tephra_grn_load"],
    "version_override": [],
    "grid_policy": "regrid",
    "regrid_n": 100,
    "domain_lon": [-100.62, -96.62],
    "domain_lat": [17.02, 21.02],
    "extra_thresholds": [],
}


# --------------------------------------------------------------------------
# Configuración
# --------------------------------------------------------------------------
def _scalar(s):
    s = s.strip()
    if s.lower() in ("true", "false"):
        return s.lower() == "true"
    for cast in (int, float):
        try:
            return cast(s)
        except ValueError:
            pass
    return s.strip("'\"")


def load_config(path):
    """Lee config.yml con PyYAML si existe; si no, con un lector mínimo; aplica valores por defecto."""
    text = Path(path).read_text(encoding="utf-8")
    try:
        import yaml  # type: ignore
        cfg = yaml.safe_load(text)
    except ImportError:
        cfg = {}
        for line in text.splitlines():
            line = line.split("#", 1)[0].rstrip()
            if not line.strip() or ":" not in line:
                continue
            key, val = line.split(":", 1)
            val = val.strip()
            if val.startswith("["):
                inner = val.strip("[]")
                cfg[key.strip()] = [_scalar(v) for v in inner.split(",") if v.strip()]
            else:
                cfg[key.strip()] = _scalar(val)
    for k, v in DEFAULTS.items():
        cfg.setdefault(k, v)
    cfg["fill_limit"] = float(cfg["fill_limit"])  # por si YAML lo lee como texto
    if cfg["grid_policy"] not in ("native", "regrid"):
        sys.exit(f"[ERROR] grid_policy debe ser 'native' o 'regrid', no {cfg['grid_policy']!r}")
    ov = {}
    for item in cfg["version_override"] or []:
        y, v = str(item).split(":")
        ov[int(y)] = v.strip()
    cfg["_version_override"] = ov
    return cfg


BASEMAP_DIR = Path(__file__).resolve().parent.parent / "data" / "basemap"


def load_basemap():
    """Polilíneas (estados, costa), etiquetas de estados y ciudades desde data/basemap/."""
    bm = {}
    for k in ("states", "coast", "state_labels", "cities"):
        f = BASEMAP_DIR / f"{k}.csv"
        bm[k] = pd.read_csv(f) if f.exists() else None
    return bm


def decorate(ax, bm, cfg, ext, labels=True):
    """Dibuja mapa base (estados CONABIO, costa Natural Earth), ciudades y cráter sobre un eje."""
    if bm is not None:
        for key, lw, col in (("states", 0.6, "0.25"), ("coast", 1.0, "k")):
            d = bm.get(key)
            if d is not None:
                for _, g in d.groupby("id"):
                    ax.plot(g["lon"], g["lat"], color=col, lw=lw, zorder=3)
        sl = bm.get("state_labels")
        if labels and sl is not None:
            for _, r in sl.iterrows():
                ax.text(r["lon"], r["lat"], r["name"], fontsize=6.5, style="italic", color="0.35",
                        ha="center", va="center", zorder=4)
        ct = bm.get("cities")
        if labels and ct is not None:
            ax.scatter(ct["lon"], ct["lat"], s=9, c="white", edgecolors="k", linewidths=0.6, zorder=5)
            for _, r in ct.iterrows():
                ax.annotate(r["name"], (r["lon"], r["lat"]), xytext=(r["dx"], r["dy"]),
                            textcoords="offset points", ha=r["ha"], fontsize=7, zorder=6,
                            bbox=dict(fc="white", ec="none", alpha=0.6, pad=0.5))
    ax.plot(cfg["crater_lon"], cfg["crater_lat"], marker="^", color="b", mec="k", ms=8, zorder=7)
    ax.set_xlim(ext[0], ext[1])
    ax.set_ylim(ext[2], ext[3])
    ax.set_aspect(1 / np.cos(np.deg2rad(19.0)))
    ax.set_xticks(np.arange(-100, -96.9, 1.0))
    ax.set_yticks(np.arange(17.5, 21.0, 0.5))


def sfx(t):
    return f"_gt{t:g}"


# --------------------------------------------------------------------------
# Procesamiento de un archivo
# --------------------------------------------------------------------------
class Accumulator:
    """Conteos por celda y métricas horarias de LOAD > umbral en una rejilla (lat, lon)."""

    def __init__(self, lat, lon, thrs):
        self.lat, self.lon, self.thrs = np.asarray(lat), np.asarray(lon), thrs
        ny, nx = len(lat), len(lon)
        self.area = cell_area_km2(self.lat, self.lon)
        self.area_flat = self.area.ravel().astype("f8")
        self.counts = {t: np.zeros((ny, nx), dtype=np.int64) for t in thrs}
        self.valid = np.zeros((ny, nx), dtype=np.int64)
        self.rec = {"max_load": [], "max_idx": []}
        for t in thrs:
            self.rec[f"nc{sfx(t)}"] = []
            self.rec[f"ar{sfx(t)}"] = []

    def update(self, a):
        n = a.shape[0]
        ok = ~np.isnan(a)
        self.valid += ok.sum(axis=0)
        for t in self.thrs:
            exc = a > np.float32(t)  # estricto, en float32 (como los datos); NaN -> False
            self.counts[t] += exc.sum(axis=0)
            flat = exc.reshape(n, -1)
            self.rec[f"nc{sfx(t)}"].append(flat.sum(axis=1))
            self.rec[f"ar{sfx(t)}"].append((flat * self.area_flat).sum(axis=1))  # sin matmul (avisos de BLAS)
        b = np.where(ok, a, -np.inf).reshape(n, -1)
        idx = b.argmax(axis=1)
        mx = b[np.arange(n), idx]
        self.rec["max_load"].append(np.where(np.isfinite(mx), mx, np.nan))
        self.rec["max_idx"].append(idx)

    def result(self, base, rec_ok):
        thr, nx = self.thrs[0], len(self.lon)
        max_idx = np.concatenate(self.rec["max_idx"])
        df = base.copy()
        df["max_load_kg_m2"] = np.concatenate(self.rec["max_load"])
        df["max_lat"] = self.lat[max_idx // nx]
        df["max_lon"] = self.lon[max_idx % nx]
        for t in self.thrs:
            df[f"ncells{'' if t == thr else sfx(t)}"] = np.concatenate(self.rec[f"nc{sfx(t)}"])
            df[f"area_km2{'' if t == thr else sfx(t)}"] = np.concatenate(self.rec[f"ar{sfx(t)}"])
        df["month"] = df["time"].dt.month
        df["hour_of_day"] = df["time"].dt.hour
        df = df[rec_ok].reset_index(drop=True)   # los registros sin ningún dato no son pronósticos
        return dict(counts=self.counts[thr], valid=self.valid, hourly=df,
                    lat=self.lat, lon=self.lon, area=self.area)


def fmt_time(ts):
    return ts.strftime("%Y-%m-%d %H:%M")


def process_file(path, scenario, year, cfg, thrs):
    nc = NCFile(path, cfg["variable_candidates"])
    version, version_full, version_src = detect_version(nc.attrs, nc.varname, cfg["_version_override"].get(year))
    hidx, tsource, n_mismatch, warns = resolve_hours(nc.horas, nc.nrec, year)
    if version == "desconocida":
        warns.append("no se pudo determinar la versión de Fall3D (usa version_override en config.yml)")
    dom_lat, dom_lon = cfg["domain_lat"], cfg["domain_lon"]
    if abs(nc.lat.min() - dom_lat[0]) > 0.01 or abs(nc.lat.max() - dom_lat[1]) > 0.01 \
            or abs(nc.lon.min() - dom_lon[0]) > 0.01 or abs(nc.lon.max() - dom_lon[1]) > 0.01:
        warns.append("el dominio del archivo difiere del configurado (domain_lat/domain_lon)")

    # Rejilla de trabajo (y nativa, si se interpola: los combinados de una sola versión se
    # calculan sin interpolar)
    rg = None
    if cfg["grid_policy"] == "regrid":
        tlat, tlon = target_axes(cfg)
        if not is_same_grid(nc.lat, nc.lon, tlat, tlon):
            rg = Regridder(nc.lat, nc.lon, tlat, tlon)
    acc_nat = Accumulator(nc.lat, nc.lon, thrs)
    acc_work = Accumulator(tlat, tlon, thrs) if rg is not None else acc_nat
    grid_native = f"{nc.nlat}x{nc.nlon}"
    grid_used = f"{len(acc_work.lat)}x{len(acc_work.lon)}"
    fill_limit, chunk = cfg["fill_limit"], int(cfg["chunk_hours"])
    rec_ok_l = []
    inv = dict(n_nan=0, n_negative=0, vmin=np.inf, vmax=-np.inf)

    for h0 in range(0, nc.nrec, chunk):
        h1 = min(h0 + chunk, nc.nrec)
        a = nc.read(h0, h1)
        n = a.shape[0]
        bad = ~np.isfinite(a) | (a >= np.float32(fill_limit))   # NaN, Inf, relleno o enmascarado
        inv["n_nan"] += int(bad.sum())
        a = np.where(bad, np.nan, a)
        ok = ~bad
        inv["n_negative"] += int((a < 0).sum())
        if ok.any():
            inv["vmin"] = min(inv["vmin"], float(np.nanmin(a)))
            inv["vmax"] = max(inv["vmax"], float(np.nanmax(a)))
        rec_ok_l.append(ok.reshape(n, -1).any(axis=1))         # registro con datos (rejilla nativa)
        acc_nat.update(a)
        if rg is not None:
            acc_work.update(rg(a))                              # interpolar LOAD antes del umbral

    nrec = nc.nrec
    fname, varname = Path(path).name, nc.varname
    nc.close()

    rec_ok = np.concatenate(rec_ok_l)
    t0 = pd.Timestamp(year, 1, 1)
    base = pd.DataFrame({"scenario": scenario, "year": year, "hour_index": hidx,
                         "time": t0 + pd.to_timedelta(hidx, unit="h")})
    work = acc_work.result(base, rec_ok)
    native = acc_nat.result(base, rec_ok) if rg is not None else work

    # Estado del año a partir de las horas con datos (no se fija a mano)
    exp = expected_hours(year)
    n_valid = int(rec_ok.sum())
    if n_valid == 0:
        status, last_hour, first_hour = "sin datos", -1, -1
    else:
        last_hour, first_hour = int(hidx[rec_ok].max()), int(hidx[rec_ok].min())
        status = "parcial" if last_hour < exp - 1 else ("con huecos" if n_valid < exp else "completo")
    n_missing = (last_hour + 1 - n_valid) if n_valid else 0
    last_time = fmt_time(t0 + pd.Timedelta(hours=last_hour)) if n_valid else ""
    first_time = fmt_time(t0 + pd.Timedelta(hours=first_hour)) if n_valid else ""
    parts = []
    if status == "parcial":
        parts.append(f"parcial (hasta {last_time})")
    if status == "sin datos":
        parts.append("sin datos")
    if n_missing > 0:
        parts.append(f"{n_missing} h sin datos")
    inv = dict(scenario=scenario, year=year, file=fname, variable=varname,
               fall3d_version=version, fall3d_version_full=version_full, version_source=version_src,
               grid_native=grid_native, grid=grid_used, regridded=rg is not None,
               n_hours=nrec, n_valid_hours=n_valid, expected_full_year=exp,
               first_time=first_time, last_time=last_time, time_source=tsource,
               n_horas_mismatch=n_mismatch, n_missing_hours=n_missing, status=status,
               complete=(status == "completo"), partial_flag=(status == "parcial"),
               nota="; ".join(parts), **inv)
    for w in warns:
        print(f"[AVISO] {fname}: {w}", file=sys.stderr)
    return dict(work=work, native=native, inv=inv)


# --------------------------------------------------------------------------
# Resúmenes
# --------------------------------------------------------------------------
def summarize_group(g, thrs):
    thr = thrs[0]
    d = {"n_hours": len(g), "max_load_kg_m2": g["max_load_kg_m2"].max()}
    for t in thrs:
        s = "" if t == thr else sfx(t)
        exc = g[f"ncells{s}"] > 0
        ar = g.loc[exc, f"area_km2{s}"]
        d[f"hours_exc{s}"] = int(exc.sum())
        d[f"frac_hours_exc{s}"] = float(exc.mean()) if len(g) else np.nan
        d[f"mean_cells_exc_hours{s}"] = g.loc[exc, f"ncells{s}"].mean() if exc.any() else 0.0
        d[f"area_km2_mean_all_hours{s}"] = g[f"area_km2{s}"].mean()
        d[f"area_km2_mean_exc_hours{s}"] = ar.mean() if exc.any() else 0.0
        d[f"area_km2_p50_exc_hours{s}"] = ar.quantile(0.5) if exc.any() else 0.0
        d[f"area_km2_p90_exc_hours{s}"] = ar.quantile(0.9) if exc.any() else 0.0
        d[f"area_km2_p99_exc_hours{s}"] = ar.quantile(0.99) if exc.any() else 0.0
        d[f"area_km2_max{s}"] = g[f"area_km2{s}"].max()
    return d


def table(H, keys, thrs):
    rows = []
    for key, g in H.groupby(keys, sort=True):
        key = key if isinstance(key, tuple) else (key,)
        row = dict(zip(keys, key))
        row.update(summarize_group(g, thrs))
        rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Periodos: años y combinados
# --------------------------------------------------------------------------
# Hay dos "espacios" de datos: "work" (rejilla de trabajo: común con grid_policy=regrid, nativa con
# native) y "native" (rejilla de cada archivo, sin interpolar). Para años no interpolados son iguales.
#   * cada año: espacio work;
#   * combinado homogéneo (una sola versión de Fall3D, p. ej. 2018-2023 con 7.2 o 2024-2025 con v8):
#     espacio native, sin interpolar (requiere que sus años compartan rejilla nativa);
#   * combinado de todos los años (mezcla de versiones): espacio work (requiere rejilla común).
def version_tag(v):
    return "v" + str(v).replace(".", "")


def build_combos(yinfo):
    years = sorted(yinfo)
    if len(years) < 2:
        return []
    vers = {y: yinfo[y]["version"] for y in years}

    def mk(kind, ys, space):
        key = "grid_native" if space == "native" else "grid"
        grids = list(dict.fromkeys(yinfo[y][key] for y in ys))
        version = "+".join(dict.fromkeys(vers[y] for y in ys))
        return dict(kind=kind, years=ys, label=f"{ys[0]}-{ys[-1]}", version=version, space=space,
                    tag=version_tag(version) if kind == "homogeneo" else "all",
                    grid=grids[0] if len(grids) == 1 else "+".join(grids),
                    grid_native="+".join(dict.fromkeys(yinfo[y]["grid_native"] for y in ys)),
                    ok=len(grids) == 1)
    if len(set(vers.values())) == 1:
        return [mk("todos", years, "work")]
    combos = []
    for v in dict.fromkeys(vers[y] for y in years):
        ys = [y for y in years if vers[y] == v]
        if len(ys) > 1:
            combos.append(mk("homogeneo", ys, "native"))
    return combos + [mk("todos", years, "work")]


def combo_note(c, scen_years, ystat):
    parts = []
    if c["kind"] == "homogeneo":
        parts.append(f"solo Fall3D {c['version']}, rejilla nativa {c['grid']}")
    elif "+" in c["version"]:
        parts.append(f"mezcla Fall3D {c['version']}: puede reflejar el cambio de versión del modelo")
    for y in c["years"]:
        if y in scen_years and ystat[y][0] != "completo":
            st, lt, nm = ystat[y]
            parts.append(f"incluye año {st} {y}" + (f" (hasta {lt})" if st == "parcial" else
                                                     f" ({nm} h sin datos)" if st == "con huecos" else ""))
    return "; ".join(parts)


def expand_periods(H, combos):
    """H: dict espacio -> DataFrame horario."""
    parts = [H["work"].assign(period=H["work"]["year"].astype(str))]
    for c in combos:
        if c["ok"]:
            d = H[c["space"]]
            parts.append(d[d["year"].isin(c["years"])].assign(period=c["label"]))
    return pd.concat(parts, ignore_index=True)


def attach_info(df, info, keys):
    out = df.merge(info, on=["scenario", "period"], how="left")
    cols = keys + ["fall3d_version", "grid", "grid_native"] + \
        [c for c in df.columns if c not in keys] + ["nota"]
    return out[cols]


# --------------------------------------------------------------------------
# Salidas
# --------------------------------------------------------------------------
def write_netcdf(path, cfg, lat, lon, Pspace, years, yinfo, sets, space_desc):
    """NetCDF de un grupo de años con la misma rejilla. sets: [(tag, años, descripción)]."""
    from scipy.io import netcdf_file
    scs = cfg["scenarios"]
    fv = np.float32(-9999.0)
    ny, nx = len(lat), len(lon)
    ce = np.zeros((len(scs), len(years), ny, nx), dtype="i4")
    ve = np.zeros_like(ce)
    for i, s in enumerate(scs):
        for j, y in enumerate(years):
            if (s, y) in Pspace["counts"]:
                ce[i, j], ve[i, j] = Pspace["counts"][(s, y)], Pspace["valid"][(s, y)]

    def prob(c, v):
        return np.where(v > 0, c / np.maximum(v, 1), fv).astype("f4")

    f = netcdf_file(str(path), "w")
    f.title = "Probabilidad de excedencia de LOAD (FALL3D, Popocatepetl)"
    f.Conventions = "CF-1.8"
    f.threshold_kg_m2 = float(cfg["threshold"])
    f.exceedance_definition = "fraction of hourly forecasts with LOAD > threshold, per cell"
    f.grid = f"{ny}x{nx}"
    f.grid_policy = cfg["grid_policy"]
    f.grid_space = space_desc
    f.years = ",".join(str(y) for y in years)
    f.fall3d_version_by_year = ",".join(f"{y}:{yinfo[y]['version']}" for y in years)
    f.partial_years = ",".join(str(y) for y in years if yinfo[y]["status"] == "parcial") or "none"
    f.years_with_gaps = ",".join(f"{y}:{yinfo[y]['n_missing']}h" for y in years if yinfo[y]["status"] == "con huecos") or "none"
    f.source = "popo-ash-exceedance"
    f.createDimension("scenario", len(scs))
    f.createDimension("year", len(years))
    f.createDimension("lat", ny)
    f.createDimension("lon", nx)

    def var(name, data, dims, typ, **attrs):
        v = f.createVariable(name, typ, dims)
        v[:] = data
        for k, val in attrs.items():
            setattr(v, k, val)
        return v

    var("scenario", np.arange(len(scs), dtype="i4"), ("scenario",), "i",
        long_name="scenario index", names=",".join(scs))
    var("column_height_km", np.array(cfg["column_height_km"], dtype="i4"), ("scenario",), "i",
        long_name="assumed eruption column height above crater", units="km")
    var("year", np.array(years, dtype="i4"), ("year",), "i", long_name="year")
    var("lat", lat.astype("f4"), ("lat",), "f", units="degrees_north", standard_name="latitude")
    var("lon", lon.astype("f4"), ("lon",), "f", units="degrees_east", standard_name="longitude")
    var("prob_exc_year", prob(ce, ve), ("scenario", "year", "lat", "lon"), "f",
        long_name="probability of exceedance per year", units="1", _FillValue=fv)
    var("count_exc_year", ce, ("scenario", "year", "lat", "lon"), "i",
        long_name="hours with LOAD > threshold", units="1")
    var("n_valid_year", ve, ("scenario", "year", "lat", "lon"), "i",
        long_name="valid (non-missing) hours", units="1")
    for tag, ys, what in sets:
        sel = [years.index(y) for y in ys]
        ca, va = ce[:, sel].sum(axis=1).astype("i4"), ve[:, sel].sum(axis=1).astype("i4")
        yl = ",".join(map(str, ys))
        var(f"prob_exc_{tag}", prob(ca, va), ("scenario", "lat", "lon"), "f",
            long_name=f"probability of exceedance, {what} ({yl})", units="1", _FillValue=fv)
        var(f"count_exc_{tag}", ca, ("scenario", "lat", "lon"), "i",
            long_name=f"hours with LOAD > threshold, {what} ({yl})", units="1")
        var(f"n_valid_{tag}", va, ("scenario", "lat", "lon"), "i", long_name=f"valid hours, {what} ({yl})", units="1")
    f.close()


def netcdf_plan(yinfo, combos):
    """Archivos NetCDF: [(nombre, espacio, rejilla, años, sets)]. El de la rejilla de trabajo lleva todos los
    años (y el combinado de todos); los años interpolados tienen además un archivo en su rejilla nativa."""
    years = sorted(yinfo)
    plan = []
    wgroups = {}
    for y in years:
        wgroups.setdefault(yinfo[y]["grid"], []).append(y)
    for g, ys in wgroups.items():
        name = "prob_exceedance.nc" if len(wgroups) == 1 else f"prob_exceedance_{g}.nc"
        plan.append([name, "work", g, ys, []])
    ngroups = {}
    for y in years:
        if yinfo[y]["regridded"]:
            ngroups.setdefault(yinfo[y]["grid_native"], []).append(y)
    for g, ys in ngroups.items():
        plan.append([f"prob_exceedance_{g}.nc", "native", g, ys, []])
    for c in combos:
        if not c["ok"]:
            continue
        what = f"all years (Fall3D {c['version']})" if c["kind"] == "todos" else f"Fall3D {c['version']} only, native grid"
        # Un combinado homogéneo sin años interpolados tiene los mismos datos en ambos espacios
        target = "work" if c["space"] == "work" or not any(yinfo[y]["regridded"] for y in c["years"]) else "native"
        for name, space, g, ys, sets in plan:
            if space == target and g == c["grid"] and set(c["years"]).issubset(ys):
                sets.append((c["tag"], c["years"], what))
                break
    return plan


def make_plots(outdir, cfg, G, P, yinfo, combos, tabs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import PowerNorm

    scs = cfg["scenarios"]
    years = sorted(yinfo)
    heights = dict(zip(scs, cfg["column_height_km"]))
    thr = cfg["threshold"]
    bm = load_basemap()
    figdir = Path(outdir) / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    foot = "Estados: CONABIO/INEGI; costa: Natural Earth"

    def ext_of(space, y):
        lat, lon = G[space][y]["lat"], G[space][y]["lon"]
        return [lon.min(), lon.max(), lat.min(), lat.max()]

    def prob(counts, valid):
        return np.where(valid > 0, counts / np.maximum(valid, 1), np.nan)

    def draw(ax, p, vmax, title, ext, cmap="YlOrRd", labels=True):
        pm = np.ma.masked_where(~(p > 0), p)
        im = ax.imshow(pm, origin="lower", extent=ext, cmap=cmap,
                       norm=PowerNorm(0.5, vmin=0, vmax=vmax), aspect="auto")
        decorate(ax, bm, cfg, ext, labels=labels)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("Lon")
        ax.set_ylabel("Lat")
        return im

    def partial_text(ys):
        txt = []
        for y in ys:
            st = yinfo[y]
            if st["status"] == "parcial":
                txt.append(f"{y} parcial (hasta {st['last_time']})")
            elif st["status"] != "completo":
                txt.append(f"{y} {st['status']}" + (f" ({st['n_missing']} h sin datos)" if st["status"] == "con huecos" else ""))
        return "; ".join(txt)

    def combo_probs(c):
        Ps = P[c["space"]]
        out = {}
        for s in scs:
            ks = [(s, y) for y in c["years"] if (s, y) in Ps["counts"]]
            if ks:
                out[s] = prob(sum(Ps["counts"][k] for k in ks), sum(Ps["valid"][k] for k in ks))
        return out

    valid_combos = [c for c in combos if c["ok"]]
    first_version = yinfo[years[0]]["version"]
    maps = valid_combos if valid_combos else [dict(kind="todos", years=[years[0]], label=str(years[0]),
                                                    version=first_version, space="work", tag="all", grid=yinfo[years[0]]["grid"])]
    has_homog = any(c["kind"] == "homogeneo" for c in maps)
    for c in maps:
        if c["kind"] == "homogeneo":
            fname = "maps_prob_combined.png" if c["version"] == first_version else f"maps_prob_combined_{c['tag']}.png"
        else:
            fname = "maps_prob_combined_all.png" if has_homog else "maps_prob_combined.png"
        pc = combo_probs(c)
        if not pc:
            continue
        ext = ext_of(c["space"], c["years"][0])
        vmax = max(np.nanmax(p) for p in pc.values())
        fig, axs = plt.subplots(1, len(pc), figsize=(6.2 * len(pc), 5.9), squeeze=False)
        for ax, s in zip(axs[0], pc):
            im = draw(ax, pc[s], vmax, f"{s} (columna {heights[s]} km), {c['label']}", ext)
        fig.colorbar(im, ax=axs.ravel().tolist(), label=f"P(LOAD > {thr:g} kg/m2)", shrink=0.85)
        sub = f"Fall3D {c['version']}, rejilla {c['grid']}" + (
            " nativa" if c["kind"] == "homogeneo" else "") + (
            " (mezcla de versiones: puede reflejar el cambio del modelo)" if "+" in c["version"] else "")
        pt = partial_text(c["years"])
        fig.suptitle(f"Probabilidad de excedencia (fracción de horas). ▲ = cráter. {foot}\n{sub}"
                     + (f". {pt}" if pt else ""), fontsize=10)
        fig.savefig(figdir / fname, dpi=150, bbox_inches="tight")
        plt.close(fig)

    Pw = P["work"]
    for s in scs:
        ys = [y for y in years if (s, y) in Pw["counts"]]
        if not ys:
            continue
        ps = {y: prob(Pw["counts"][(s, y)], Pw["valid"][(s, y)]) for y in ys}
        vm = max(np.nanmax(p) for p in ps.values())
        ncol = 4
        nrow = int(np.ceil(len(ys) / ncol))
        fig, axs = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.8 * nrow), squeeze=False)
        for ax in axs.ravel():
            ax.axis("off")
        for ax, y in zip(axs.ravel(), ys):
            ax.axis("on")
            st = yinfo[y]
            extra = (f"\nparcial, hasta {st['last_time']}" if st["status"] == "parcial" else
                     f"\ncon huecos ({st['n_missing']} h sin datos)" if st["status"] == "con huecos" else
                     f"\n{st['status']}" if st["status"] != "completo" else "")
            interp = ", interpolado" if st["regridded"] else ""
            im = draw(ax, ps[y], vm, f"{s} {y} (Fall3D {st['version']}{interp}){extra}", ext_of("work", y), labels=False)
        fig.colorbar(im, ax=axs.ravel().tolist(), label=f"P(LOAD > {thr:g} kg/m2)", shrink=0.8)
        fig.savefig(figdir / f"maps_prob_by_year_{s}.png", dpi=130, bbox_inches="tight")
        plt.close(fig)

    c0 = maps[0]
    pc = combo_probs(c0)
    pairs = [(a, b) for a, b in [("c05M", "c03M"), ("c10M", "c05M"), ("c10M", "c03M")] if a in pc and b in pc]
    if pairs:
        ext = ext_of(c0["space"], c0["years"][0])
        fig, axs = plt.subplots(1, len(pairs), figsize=(6.2 * len(pairs), 5.9), squeeze=False)
        dm = max(np.nanmax(np.abs(pc[a] - pc[b])) for a, b in pairs)
        for ax, (a, b) in zip(axs[0], pairs):
            im = ax.imshow(pc[a] - pc[b], origin="lower", extent=ext, cmap="RdBu_r",
                           vmin=-dm, vmax=dm, aspect="auto")
            decorate(ax, bm, cfg, ext)
            ax.set_title(f"{a} - {b}", fontsize=9)
            ax.set_xlabel("Lon")
            ax.set_ylabel("Lat")
        fig.colorbar(im, ax=axs.ravel().tolist(), label="Diferencia de probabilidad", shrink=0.85)
        fig.suptitle(f"Diferencia entre escenarios, {c0['label']} (Fall3D {c0['version']}). {foot}", fontsize=10)
        fig.savefig(figdir / "maps_prob_difference.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    # Barras: fracción de horas con excedencia por año
    ty = tabs["by_year"]
    ty = ty[ty["period"].isin([str(y) for y in years])]
    fig, ax = plt.subplots(figsize=(max(9, 1.2 * len(years) + 3), 4.8))
    w = 0.8 / len(scs)
    for i, s in enumerate(scs):
        d = ty[ty["scenario"] == s]
        xs = [years.index(int(p)) + i * w for p in d["period"]]
        bars = ax.bar(xs, 100 * d["frac_hours_exc"], w, label=f"{s} ({heights[s]} km)")
        for bar, p in zip(bars, d["period"]):
            if yinfo[int(p)]["status"] != "completo":
                bar.set_hatch("//")
    ax.set_xticks([j + 0.4 - w / 2 for j in range(len(years))])
    ax.set_xticklabels([f"{y}{'*' if yinfo[y]['status'] != 'completo' else ''}\nFall3D {yinfo[y]['version']}" for y in years],
                       fontsize=8)
    ax.set_ylabel(f"% de horas con LOAD > {thr:g} kg/m2 en alguna celda")
    ax.set_title("Horas con excedencia por año (* = parcial o con huecos, achurado)")
    pt = partial_text(years)
    if pt:
        fig.text(0.01, -0.02, pt, fontsize=8, ha="left", va="top")
    ax.legend()
    fig.savefig(figdir / "hours_exceedance_by_year.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    H = tabs["hourly"]
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    data, labs = [], []
    for s in scs:
        v = H.loc[(H["scenario"] == s) & (H["ncells"] > 0), "area_km2"].values
        if len(v):
            data.append(v)
            labs.append(f"{s}\n(n={len(v)})")
    if data:
        ax.boxplot(data, tick_labels=labs, showfliers=True)
        ax.set_yscale("log")
        ax.set_ylabel("Área con LOAD > umbral (km2), horas con excedencia")
        ax.set_title("Área de influencia por escenario (todos los años, rejilla de trabajo)")
        fig.savefig(figdir / "area_by_scenario.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    for key, xcol, fname, xl in [("monthly", "month", "monthly_exceedance.png", "Mes"),
                                  ("diurnal", "hour_of_day", "diurnal_exceedance.png", "Hora del día")]:
        d = tabs[key]
        fig, ax = plt.subplots(figsize=(7.5, 4.2))
        for s in scs:
            dd = d[d["scenario"] == s].sort_values(xcol)
            ax.plot(dd[xcol], 100 * dd["frac_hours_exc"], marker="o", label=s)
        ax.set_xlabel(xl)
        ax.set_ylabel("% de horas con excedencia")
        ax.set_title(f"{c0['label']} (Fall3D {c0['version']})", fontsize=9)
        ax.legend()
        ax.grid(alpha=0.3)
        fig.savefig(figdir / fname, dpi=150, bbox_inches="tight")
        plt.close(fig)


# --------------------------------------------------------------------------
# Principal
# --------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--config", default="config/config.yml")
    ap.add_argument("--data-dir")
    ap.add_argument("--output-dir")
    ap.add_argument("--scenarios", nargs="+")
    ap.add_argument("--years", nargs="+", type=int)
    ap.add_argument("--threshold", type=float)
    ap.add_argument("--grid-policy", choices=["native", "regrid"])
    ap.add_argument("--skip-plots", action="store_true")
    a = ap.parse_args(argv)

    cfg = load_config(a.config)
    # Prioridad de data_dir: --data-dir > variable de entorno POPO_DATA_DIR > config.yml
    if a.data_dir:
        cfg["data_dir"] = a.data_dir
    elif os.environ.get("POPO_DATA_DIR"):
        cfg["data_dir"] = os.environ["POPO_DATA_DIR"]
    if a.output_dir:
        cfg["output_dir"] = a.output_dir
    if a.threshold is not None:
        cfg["threshold"] = a.threshold
    if a.grid_policy:
        cfg["grid_policy"] = a.grid_policy
    if a.scenarios:
        keep = [s for s in cfg["scenarios"] if s in a.scenarios]
        cfg["column_height_km"] = [h for s, h in zip(cfg["scenarios"], cfg["column_height_km"]) if s in keep]
        cfg["scenarios"] = keep
    if a.years:
        cfg["years"] = sorted(a.years)
    thrs = [float(cfg["threshold"])] + [float(t) for t in cfg["extra_thresholds"]]

    if not Path(str(cfg["data_dir"])).is_dir():
        sys.exit(f"[ERROR] data_dir no existe: {cfg['data_dir']!r}. Edita config/config.yml, "
                 "usa --data-dir o define la variable de entorno POPO_DATA_DIR.")
    outdir = Path(cfg["output_dir"])
    (outdir / "tables").mkdir(parents=True, exist_ok=True)
    (outdir / "netcdf").mkdir(parents=True, exist_ok=True)

    SP = ("work", "native")
    P = {sp: dict(counts={}, valid={}) for sp in SP}
    G = {sp: {} for sp in SP}
    hourly = {sp: [] for sp in SP}
    inv = []
    t0 = time.time()
    for s in cfg["scenarios"]:
        for y in cfg["years"]:
            path = Path(cfg["data_dir"]) / cfg["file_pattern"].format(scenario=s, year=y)
            if not path.exists():
                print(f"[AVISO] no existe {path.name}; se omite", file=sys.stderr)
                continue
            r = process_file(path, s, y, cfg, thrs)
            for sp in SP:
                P[sp]["counts"][(s, y)] = r[sp]["counts"]
                P[sp]["valid"][(s, y)] = r[sp]["valid"]
                hourly[sp].append(r[sp]["hourly"])
                G[sp].setdefault(y, dict(lat=r[sp]["lat"], lon=r[sp]["lon"], area=r[sp]["area"]))
            inv.append(r["inv"])
            i = r["inv"]
            print(f"[{time.time() - t0:6.1f}s] {path.name}: Fall3D {i['fall3d_version']}, {i['variable']}, "
                  f"{i['grid_native']}->{i['grid']}, {i['n_hours']} registros, {i['status']}, "
                  f"horas con exc={(r['work']['hourly']['ncells'] > 0).sum()}, max={i['vmax']:.3f} kg/m2", flush=True)
    if not inv:
        sys.exit("No se procesó ningún archivo.")

    inv_df = pd.DataFrame(inv)
    # Información por año (se verifica que los escenarios coincidan en versión y rejilla)
    yinfo = {}
    for y, g in inv_df.groupby("year"):
        first = g.iloc[0]
        for col in ("fall3d_version", "grid", "grid_native"):
            if g[col].nunique() > 1:
                print(f"[AVISO] {y}: los escenarios difieren en {col}: {g[col].tolist()}", file=sys.stderr)
        worst = g.iloc[(g["status"] != "completo").argmax()] if (g["status"] != "completo").any() else first
        yinfo[int(y)] = dict(version=first["fall3d_version"], grid=first["grid"], grid_native=first["grid_native"],
                             regridded=bool(first["regridded"]), status=worst["status"],
                             last_time=worst["last_time"], n_missing=int(worst["n_missing_hours"]))
    combos = build_combos(yinfo)
    for c in combos:
        if not c["ok"]:
            print(f"[AVISO] periodo {c['label']} ({c['kind']}) omitido: sus años usan rejillas distintas "
                  f"({c['grid']}); usa grid_policy: regrid", file=sys.stderr)

    # Información (versión, rejilla, nota) por escenario y periodo
    rows = []
    for _, r in inv_df.iterrows():
        rows.append(dict(scenario=r["scenario"], period=str(r["year"]), fall3d_version=r["fall3d_version"],
                         grid=r["grid"], grid_native=r["grid_native"], nota=r["nota"]))
    ystat = {s: {int(r["year"]): (r["status"], r["last_time"], int(r["n_missing_hours"]))
                 for _, r in inv_df[inv_df.scenario == s].iterrows()} for s in cfg["scenarios"]}
    for c in combos:
        if not c["ok"]:
            continue
        for s in cfg["scenarios"]:
            if ystat.get(s):
                rows.append(dict(scenario=s, period=c["label"], fall3d_version=c["version"], grid=c["grid"],
                                 grid_native=c["grid_native"], nota=combo_note(c, set(ystat[s]), ystat[s])))
    info = pd.DataFrame(rows)

    H = {sp: pd.concat(hourly[sp], ignore_index=True) for sp in SP}
    Hp = expand_periods(H, combos)
    by_year = attach_info(table(Hp, ["scenario", "period"], thrs), info, ["scenario", "period"])
    monthly = attach_info(table(Hp, ["scenario", "period", "month"], thrs), info, ["scenario", "period", "month"])
    wet = set(cfg["wet_months"])
    Hs = Hp.assign(season=np.where(Hp["month"].isin(wet), "lluvias(may-oct)", "secas(nov-abr)"))
    seasonal = attach_info(table(Hs, ["scenario", "period", "season"], thrs), info, ["scenario", "period", "season"])
    diurnal = attach_info(table(Hp, ["scenario", "period", "hour_of_day"], thrs), info, ["scenario", "period", "hour_of_day"])

    # Probabilidad máxima y área por nivel de probabilidad
    lvl_rows, max_rows = [], []
    for s in cfg["scenarios"]:
        pers = {}
        for y in sorted(G["work"]):
            if (s, y) in P["work"]["counts"]:
                pers[str(y)] = (P["work"]["counts"][(s, y)], P["work"]["valid"][(s, y)], "work", y)
        for c in combos:
            Ps = P[c["space"]]
            ys = [y for y in c["years"] if (s, y) in Ps["counts"]]
            if c["ok"] and ys:
                pers[c["label"]] = (sum(Ps["counts"][(s, y)] for y in ys), sum(Ps["valid"][(s, y)] for y in ys),
                                    c["space"], ys[0])
        for p, (cnt, val, sp, y0) in pers.items():
            lat, lon, area = G[sp][y0]["lat"], G[sp][y0]["lon"], G[sp][y0]["area"]
            pr = np.where(val > 0, cnt / np.maximum(val, 1), np.nan)
            k = int(np.nanargmax(pr))
            max_rows.append(dict(scenario=s, period=p, max_prob=float(np.nanmax(pr)),
                                 max_prob_lat=float(lat[k // len(lon)]), max_prob_lon=float(lon[k % len(lon)]),
                                 cells_ever_exceeded=int((cnt > 0).sum()),
                                 area_km2_ever_exceeded=float(area[cnt > 0].sum())))
            for lv in cfg["prob_levels"]:
                m = np.nan_to_num(pr, nan=0.0) >= lv
                lvl_rows.append(dict(scenario=s, period=p, prob_level=lv,
                                     n_cells=int(m.sum()), area_km2=float(area[m].sum())))
    max_prob = attach_info(pd.DataFrame(max_rows), info, ["scenario", "period"])
    area_levels = attach_info(pd.DataFrame(lvl_rows), info, ["scenario", "period", "prob_level"])

    Hw = H["work"]
    exc_hours = Hw.loc[Hw["ncells"] > 0, ["scenario", "year", "hour_index", "time", "ncells",
                                            "area_km2", "max_load_kg_m2", "max_lat", "max_lon"]]
    pdefs = pd.DataFrame([dict(period=str(y), kind="año", years=str(y), fall3d_version=yinfo[y]["version"],
                               grid=yinfo[y]["grid"], grid_native=yinfo[y]["grid_native"], status=yinfo[y]["status"])
                          for y in sorted(yinfo)] +
                         [dict(period=c["label"], kind=c["kind"], years=" ".join(map(str, c["years"])),
                               fall3d_version=c["version"], grid=c["grid"], grid_native=c["grid_native"],
                               status="combinado" if c["ok"] else "omitido (rejillas distintas)") for c in combos])

    tdir = outdir / "tables"
    inv_df.to_csv(tdir / "file_inventory.csv", index=False)
    pdefs.to_csv(tdir / "period_definitions.csv", index=False)
    by_year.to_csv(tdir / "summary_scenario_period.csv", index=False)
    monthly.to_csv(tdir / "summary_monthly.csv", index=False)
    seasonal.to_csv(tdir / "summary_seasonal.csv", index=False)
    diurnal.to_csv(tdir / "summary_diurnal.csv", index=False)
    max_prob.to_csv(tdir / "max_probability.csv", index=False)
    area_levels.to_csv(tdir / "area_by_probability_level.csv", index=False)
    exc_hours.to_csv(tdir / "exceedance_hours.csv", index=False)
    Hw.to_csv(tdir / "hourly_metrics.csv.gz", index=False)

    for name, sp, g, ys, sets in netcdf_plan(yinfo, combos):
        desc = "native grid (not interpolated)" if sp == "native" else \
            ("common grid (regrid)" if cfg["grid_policy"] == "regrid" else "native grid")
        write_netcdf(outdir / "netcdf" / name, cfg, G[sp][ys[0]]["lat"], G[sp][ys[0]]["lon"], P[sp], ys, yinfo, sets, desc)

    if not a.skip_plots:
        valid_combos = [c for c in combos if c["ok"]]
        primary = valid_combos[0]["label"] if valid_combos else str(sorted(yinfo)[0])
        make_plots(outdir, cfg, G, P, yinfo, combos,
                   dict(by_year=by_year, monthly=monthly[monthly["period"] == primary],
                        diurnal=diurnal[diurnal["period"] == primary], hourly=Hw))
    print(f"Listo en {time.time() - t0:.0f} s. Salidas en {outdir.resolve()}")


if __name__ == "__main__":
    main()
