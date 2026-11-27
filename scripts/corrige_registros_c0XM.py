#!/usr/bin/env python3
"""Reemplaza registros de un c<CASO>-<ANIO>.nc (Fall3D 7.1 con TITLE "7.2", NetCDF clasico) por LOAD nuevo, sin tocar el original.

Copia el original a un archivo NUEVO y reescribe solo los registros indicados (LOAD(hora,time,lat,lon); registro i = AAAA-01-01 00:00 + i horas).
No cambia atributos, formato, dimensiones ni otras variables. Verifica al final bit a bit: registros sin cambio identicos al original,
registros reemplazados identicos al LOAD dado, y que los bytes distintos caen solo en los registros reemplazados.

El archivo de reemplazos es un .npz con:  registros (int, n)  y  LOAD (float32, n x ny x nx)  = LOAD del ULTIMO paso de time
de la simulacion correspondiente (equivale a  ncks -d time,-1 -v LOAD,lat,lon  de crea_IA.sh).

Uso:
    python scripts/corrige_registros_c0XM.py --original c03M-2018.nc --reemplazos rep_c03M-2018.npz --salida c03M-2018.corregido.nc
Los datos (.nc, .npz) no se versionan.
"""
import argparse, shutil, sys
from pathlib import Path
import numpy as np
import netCDF4 as nc


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--original", required=True); ap.add_argument("--reemplazos", required=True); ap.add_argument("--salida", required=True)
    a = ap.parse_args()
    src, dst = Path(a.original), Path(a.salida)
    if dst.exists() or dst.resolve() == src.resolve():
        sys.exit(f"[ERROR] la salida debe ser un archivo nuevo: {dst}")
    z = np.load(a.reemplazos); reg = np.asarray(z["registros"], dtype=int); new = np.asarray(z["LOAD"], dtype=np.float32)
    if len(set(reg.tolist())) != len(reg) or new.shape[0] != len(reg):
        sys.exit("[ERROR] registros repetidos o numero de campos distinto del de registros")
    with nc.Dataset(src) as o:
        if o.data_model not in ("NETCDF3_CLASSIC", "NETCDF3_64BIT_OFFSET"): sys.exit(f"[ERROR] formato inesperado: {o.data_model}")
        nrec, shp = o["LOAD"].shape[0], o["LOAD"].shape[2:]
    if reg.min() < 0 or reg.max() >= nrec or new.shape[1:] != tuple(shp):
        sys.exit(f"[ERROR] registros fuera de rango o rejilla distinta ({new.shape[1:]} vs {shp})")
    shutil.copyfile(src, dst)
    with nc.Dataset(dst, "r+") as d:
        d.set_auto_mask(False)
        for i, f in zip(reg, new): d["LOAD"][int(i), 0, :, :] = f
    # verificacion bit a bit
    rep = {int(i): k for k, i in enumerate(reg)}; bad = 0
    with nc.Dataset(src) as o, nc.Dataset(dst) as n:
        o.set_auto_mask(False); n.set_auto_mask(False)
        assert o.data_model == n.data_model and list(o.variables) == list(n.variables)
        assert all(np.array_equal(o[v][:], n[v][:]) for v in o.variables if v != "LOAD"), "variables distintas de LOAD cambiaron"
        assert {k: str(o.getncattr(k)) for k in o.ncattrs()} == {k: str(n.getncattr(k)) for k in n.ncattrs()}, "atributos globales cambiaron"
        for s in range(0, nrec, 500):
            ao, an = o["LOAD"][s:s + 500].view(np.uint32), n["LOAD"][s:s + 500].view(np.uint32)
            for j in range(ao.shape[0]):
                i = s + j
                exp = new[rep[i]].view(np.uint32).reshape(an[j].shape) if i in rep else ao[j]
                bad += not np.array_equal(an[j], exp)
    if bad: sys.exit(f"[ERROR] {bad} registros no coinciden con lo esperado; borre {dst}")
    print(f"OK: {dst} ({len(reg)} registros reemplazados de {nrec}); verificado bit a bit.")


if __name__ == "__main__":
    main()
