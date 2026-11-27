#!/usr/bin/env python3
"""Compara las salidas de Python y R (tablas CSV y NetCDF) y revisa expectativas del test sintético.

Autor: Dr. José Agustín García Reynoso (ICAyCC, UNAM) <agustin@atmosfera.unam.mx>. Licencia: MIT.
Uso: python tests/compare_outputs.py DIR_PYTHON DIR_R [--expect dataset1_regrid|dataset1_native|dataset2]
Sale con código 1 si alguna comparación o expectativa falla.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import netcdf_file

TOL = 1e-6          # tolerancia relativa para valores reales (áreas, máximos)
KEYS = {
    "file_inventory.csv": ["scenario", "year"],
    "period_definitions.csv": ["period"],
    "summary_scenario_period.csv": ["scenario", "period"],
    "summary_monthly.csv": ["scenario", "period", "month"],
    "summary_seasonal.csv": ["scenario", "period", "season"],
    "summary_diurnal.csv": ["scenario", "period", "hour_of_day"],
    "max_probability.csv": ["scenario", "period"],
    "area_by_probability_level.csv": ["scenario", "period", "prob_level"],
    "exceedance_hours.csv": ["scenario", "year", "hour_index"],
    "hourly_metrics.csv.gz": ["scenario", "year", "hour_index"],
}
fails = []


def check(cond, msg):
    print(("  OK   " if cond else "  FALLA ") + msg)
    if not cond:
        fails.append(msg)


def read(path, keys):
    d = pd.read_csv(path, keep_default_na=False, na_values=["NA", "NaN", "nan"])
    for k in keys:
        d[k] = d[k].astype(str)
    return d.sort_values(keys).reset_index(drop=True)


def norm(s):
    s = s.astype(str).str.strip()
    return s.replace({"True": "TRUE", "False": "FALSE", "nan": "", "NA": ""})


def compare_tables(dp, dr):
    print("Tablas CSV (Python vs R):")
    for name, keys in KEYS.items():
        a, b = read(dp / "tables" / name, keys), read(dr / "tables" / name, keys)
        if list(a.columns) != list(b.columns):
            check(False, f"{name}: columnas distintas {set(a.columns) ^ set(b.columns)}")
            continue
        if len(a) != len(b) or not (a[keys].values == b[keys].values).all():
            check(False, f"{name}: filas distintas ({len(a)} vs {len(b)})")
            continue
        worst, worst_col, txt_bad = 0.0, "", []
        for c in a.columns:
            if c in keys:
                continue
            is_bool = a[c].dtype == bool or b[c].dtype == bool
            x = pd.to_numeric(a[c], errors="coerce") if not is_bool else pd.Series(dtype="f8")
            y = pd.to_numeric(b[c], errors="coerce") if not is_bool else pd.Series(dtype="f8")
            if not is_bool and x.notna().all() and y.notna().all() and a[c].astype(str).ne("").all():
                x, y = x.astype("f8"), y.astype("f8")
                rel = (np.abs(x - y) / np.maximum(np.abs(x), 1e-9)).max() if len(x) else 0.0
                if rel > worst:
                    worst, worst_col = float(rel), c
            elif not (norm(a[c]) == norm(b[c])).all():
                txt_bad.append(c)
        check(worst <= TOL and not txt_bad,
              f"{name}: {len(a)} filas; dif. relativa máx {worst:.2e} ({worst_col or '-'})"
              + (f"; texto distinto en {txt_bad}" if txt_bad else ""))


def compare_netcdf(dp, dr):
    print("NetCDF (Python vs R):")
    fp = sorted(p.name for p in (dp / "netcdf").glob("*.nc"))
    fr = sorted(p.name for p in (dr / "netcdf").glob("*.nc"))
    check(fp == fr, f"mismos archivos: {fp}")
    for name in fp:
        if name not in fr:
            continue
        a, b = netcdf_file(dp / "netcdf" / name, mmap=False), netcdf_file(dr / "netcdf" / name, mmap=False)
        vs = [v for v in a.variables if v.startswith(("prob_exc", "count_exc", "n_valid"))]
        check(set(vs) == {v for v in b.variables if v.startswith(("prob_exc", "count_exc", "n_valid"))},
              f"{name}: mismas variables {sorted(vs)}")
        for v in vs:
            x = np.array(a.variables[v][:], dtype="f8")
            y = np.array(b.variables[v][:], dtype="f8")
            # Python guarda (scenario, year, lat, lon); ncdf4 en R escribe el mismo orden en disco
            check(x.shape == y.shape and np.abs(x - y).max() < 1e-6, f"{name}:{v} {x.shape} dif. máx "
                  f"{np.abs(x - y).max() if x.shape == y.shape else 'forma distinta'}")
        for att in ("grid", "years", "fall3d_version_by_year", "partial_years", "years_with_gaps", "grid_space"):
            pa, ra = getattr(a, att, b"?"), getattr(b, att, b"?")
            pa = pa.decode() if isinstance(pa, bytes) else str(pa)
            ra = ra.decode() if isinstance(ra, bytes) else str(ra)
            check(pa == ra, f"{name}: atributo {att} = {pa!r} (R {ra!r})")


def expectations(dp, case):
    print(f"Expectativas del caso {case} (Python):")
    inv = pd.read_csv(dp / "tables" / "file_inventory.csv", keep_default_na=False,
                      dtype={"fall3d_version": str, "fall3d_version_full": str}).set_index("year")
    per = pd.read_csv(dp / "tables" / "period_definitions.csv", keep_default_na=False, dtype=str).set_index("period")
    summ = pd.read_csv(dp / "tables" / "summary_scenario_period.csv", keep_default_na=False)
    eh = pd.read_csv(dp / "tables" / "exceedance_hours.csv")
    nc = sorted(p.name for p in (dp / "netcdf").glob("*.nc"))
    if case.startswith("dataset1"):
        check(sorted(inv.index) == [2021, 2022, 2024, 2025], "2023 omitido y .prev_20240424 no recogido")
        check(inv.loc[2024, "variable"] == "tephra_grn_load" and inv.loc[2025, "variable"] == "LOAD",
              "variable detectada: 2024 tephra_grn_load, 2025 LOAD")
        check(inv.loc[2021, "fall3d_version"] == "7.2" and str(inv.loc[2024, "fall3d_version"]) == "8"
              and str(inv.loc[2025, "fall3d_version_full"]) == "8.2.1"
              and str(inv.loc[2024, "fall3d_version_full"]) == "8.0.1", "versiones 7.2 / 8.0.1 / 8.2.1")
        check(inv.loc[2021, "time_source"] == "index" and inv.loc[2024, "time_source"] == "horas",
              "tiempo: 2021 por índice ('horas' escalar 0), 2024 por 'horas'")
        check(int(inv.loc[2024, "n_horas_mismatch"]) == 30 and int(inv.loc[2024, "n_valid_hours"]) == 69
              and int(inv.loc[2024, "n_missing_hours"]) == 5, "2024: hueco de 4 h + 1 registro relleno = 5 h sin datos")
        check(inv.loc[2024, "last_time"] == "2024-01-04 01:00", f"2024: última hora real {inv.loc[2024, 'last_time']}")
        h24 = eh[eh.year == 2024]
        check(h24.hour_index.max() == 73 and not h24.hour_index.between(40, 43).any(),
              "exceedance_hours 2024 usa la hora real (sin 40-43, hasta 73)")
        check((inv["status"] == "parcial").all(), "todos los años marcados parciales (72 h)")
        if case == "dataset1_regrid":
            check((inv["grid"] == "100x100").all() and bool(inv.loc[2024, "regridded"]) and not bool(inv.loc[2021, "regridded"]),
                  "regrid: todo a 100x100; solo v8 interpolado")
            check(set(per.index) == {"2021", "2022", "2024", "2025", "2021-2022", "2024-2025", "2021-2025"}
                  and per.loc["2021-2022", "kind"] == "homogeneo" and per.loc["2021-2025", "fall3d_version"] == "7.2+8",
                  "periodos: 2021-2022 homogéneo (7.2), 2024-2025 homogéneo (v8) y 2021-2025 todos (7.2+8)")
            check(per.loc["2024-2025", "grid"] == "101x101" and per.loc["2024", "grid"] == "100x100",
                  "regrid: años v8 en 100x100, combinado v8 en su rejilla nativa 101x101")
            v8_native_combo(dp, "prob_exceedance_101x101.nc", "prob_exceedance.nc", 100)
            nota = summ.set_index("period").loc["2021-2025", "nota"]
            check("cambio de versión" in nota, f"nota del combinado mixto: {nota!r}")
            check(nc == ["prob_exceedance.nc", "prob_exceedance_101x101.nc"], f"NetCDF común + nativo v8: {nc}")
        else:
            check(inv.loc[2021, "grid"] == "100x100" and inv.loc[2024, "grid"] == "101x101", "native: cada año en su rejilla")
            check(per.loc["2021-2025", "status"].startswith("omitido") and per.loc["2021-2022", "status"] == "combinado"
                  and per.loc["2024-2025", "status"] == "combinado",
                  "native: combinado de todos omitido; 2021-2022 y 2024-2025 combinados")
            v8_native_combo(dp, "prob_exceedance_101x101.nc", "prob_exceedance_100x100.nc", None)
            check(nc == ["prob_exceedance_100x100.nc", "prob_exceedance_101x101.nc"], f"NetCDF por rejilla: {nc}")
    if case == "dataset2":
        check(inv.loc[2023, "status"] == "completo" and inv.loc[2024, "status"] == "completo"
              and inv.loc[2025, "status"] == "parcial", "2023 y 2024 (8784 h) completos; 2025 parcial")
        check(inv.loc[2024, "nota"] == "" and "hasta 2025-07-28 07:00" in inv.loc[2025, "nota"],
              f"notas: 2024 {inv.loc[2024, 'nota']!r}; 2025 {inv.loc[2025, 'nota']!r}")
        check(per.loc["2023-2025", "fall3d_version"] == "7.2+8" and "2023-2023" not in per.index
              and per.loc["2024-2025", "grid"] == "11x11",
              "con un solo año 7.2 no hay combinado 7.2; sí 2024-2025 v8 (11x11 nativo) y 2023-2025 mixto")
        check(nc == ["prob_exceedance.nc", "prob_exceedance_11x11.nc"], f"NetCDF común + nativo v8: {nc}")
        v8_native_combo(dp, "prob_exceedance_11x11.nc", "prob_exceedance.nc", 10)


def v8_native_combo(dp, native_file, work_file, work_n):
    """El combinado v8 está en la rejilla nativa y es la suma de los conteos nativos por año."""
    a = netcdf_file(dp / "netcdf" / native_file, mmap=False)
    vs = set(a.variables)
    check("prob_exc_v8" in vs and "prob_exc_all" not in vs, f"{native_file}: tiene prob_exc_v8 y no prob_exc_all")
    ce = np.array(a.variables["count_exc_year"][:], "f8")
    cv = np.array(a.variables["count_exc_v8"][:], "f8")
    gs = getattr(a, "grid_space", b"")
    gs = gs.decode() if isinstance(gs, bytes) else str(gs)
    check(np.array_equal(ce.sum(axis=1), cv) and "native" in gs,
          f"{native_file}: count_exc_v8 = suma de conteos nativos por año; grid_space={gs!r}")
    w = netcdf_file(dp / "netcdf" / work_file, mmap=False)
    check("prob_exc_v8" not in w.variables, f"{work_file}: sin prob_exc_v8 (va en la rejilla nativa)")
    if work_n:
        check("prob_exc_all" in w.variables and np.array(w.variables["lat"][:]).size == work_n,
              f"{work_file}: rejilla común {work_n}x{work_n} con prob_exc_all")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir_python")
    ap.add_argument("dir_r")
    ap.add_argument("--expect")
    a = ap.parse_args()
    dp, dr = Path(a.dir_python), Path(a.dir_r)
    compare_tables(dp, dr)
    compare_netcdf(dp, dr)
    if a.expect:
        expectations(dp, a.expect)
    print(f"Resultado: {'FALLAS: ' + str(len(fails)) if fails else 'todo coincide'}")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
