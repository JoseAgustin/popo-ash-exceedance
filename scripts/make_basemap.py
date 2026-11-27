#!/usr/bin/env python3
"""Genera data/basemap/*.csv (polilíneas recortadas al dominio) desde Natural Earth.

Los GeoJSON originales (1:10m, dominio público) se descargan a data_cache/ (ignorado por git):
  python scripts/make_basemap.py --download
Luego se recortan al dominio (con margen) y se guardan como CSV pequeños (id, lon, lat).

Los estados pueden venir de un shapefile propio (polígonos, lon/lat WGS84) en lugar de Natural Earth:
  python scripts/make_basemap.py --states-shp RUTA/destdv1gw.shp --name-field ENTIDAD
En ese caso solo se versionan los CSV recortados; el shapefile completo no se copia al repositorio.
La costa sigue viniendo de Natural Earth. Solo requiere la biblioteca estándar de Python.

Autor: Dr. José Agustín García Reynoso (ICAyCC, UNAM) <agustin@atmosfera.unam.mx>
Licencia: MIT (ver LICENSE). Los datos de estados (CONABIO) conservan sus propias condiciones.
"""
import argparse
import csv
import json
import struct
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data_cache"
OUT = ROOT / "data" / "basemap"
BASE_URL = "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/"
LAYERS = {
    "states": "ne_10m_admin_1_states_provinces",
    "coast": "ne_10m_coastline",
}
DOM = (-100.62, -96.62, 17.02, 21.02)  # dominio de los datos
# Estados sin etiqueta propia (muy pequeños o con ciudad del mismo nombre) y renombres
SKIP_LABELS = {"Distrito Federal", "Tlaxcala"}
RENAME = {"México": "Edo. de México"}  # el nombre mostrado ya es "México" en ambas fuentes
# Posiciones manuales de etiquetas para evitar choques con ciudades (lon, lat)
LABEL_POS = {"México": (-99.85, 19.8), "Puebla": (-97.95, 18.45), "Morelos": (-99.05, 18.55),
             "Guanajuato": (-100.2, 20.75), "Querétaro": (-99.8, 20.35),
             "Michoacán": (-100.2, 19.0), "Veracruz": (-97.65, 20.35)}
BOX = (-101.3, -95.9, 16.4, 21.7)  # lon_min, lon_max, lat_min, lat_max (dominio + margen)


# Nombres de estado para mostrar (shapefile CONABIO usa mayúsculas sin acentos)
DISPLAY_NAMES = {
    "MEXICO": "México", "MICHOACAN DE OCAMPO": "Michoacán", "QUERETARO DE ARTEAGA": "Querétaro",
    "VERACRUZ DE IGNACIO DE LA LLAVE": "Veracruz", "SAN LUIS POTOSI": "San Luis Potosí",
    "COAHUILA DE ZARAGOZA": "Coahuila", "NUEVO LEON": "Nuevo León", "YUCATAN": "Yucatán",
}


def read_dbf(path, encoding="cp850"):
    """Lee un .dbf (dBase III) con la biblioteca estándar. Devuelve lista de dicts."""
    b = Path(path).read_bytes()
    nrec, hlen, rlen = struct.unpack("<IHH", b[4:12])
    fields, i = [], 32
    while b[i] != 0x0D:
        fields.append((b[i:i + 11].split(b"\0")[0].decode("ascii"), b[i + 16]))
        i += 32
    rows = []
    for r in range(nrec):
        rec = b[hlen + r * rlen: hlen + (r + 1) * rlen]
        o, row = 1, {}
        for name, ln in fields:
            row[name] = rec[o:o + ln].decode(encoding).strip()
            o += ln
        rows.append(row)
    return rows


def read_shp_polygons(path):
    """Lee un .shp de polígonos (tipo 5). Devuelve lista (una por registro) de anillos [(x, y), ...]."""
    b = Path(path).read_bytes()
    if struct.unpack("<i", b[32:36])[0] != 5:
        raise ValueError("Se esperaba un shapefile de polígonos (tipo 5)")
    pos, out = 100, []
    while pos < len(b):
        clen = struct.unpack(">i", b[pos + 4:pos + 8])[0] * 2
        c = b[pos + 8:pos + 8 + clen]
        pos += 8 + clen
        if struct.unpack("<i", c[:4])[0] == 0:  # registro nulo
            out.append([])
            continue
        nparts, npts = struct.unpack("<ii", c[36:44])
        parts = list(struct.unpack(f"<{nparts}i", c[44:44 + 4 * nparts])) + [npts]
        o = 44 + 4 * nparts
        pts = [struct.unpack("<2d", c[o + 16 * k:o + 16 * k + 16]) for k in range(npts)]
        out.append([pts[parts[k]:parts[k + 1]] for k in range(nparts)])
    return out


def lines_from_geometry(geom):
    t, c = geom["type"], geom["coordinates"]
    if t == "LineString":
        return [c]
    if t == "MultiLineString":
        return list(c)
    if t == "Polygon":
        return list(c)
    if t == "MultiPolygon":
        return [r for poly in c for r in poly]
    return []


def clip_runs(line):
    """Tramos de la polilínea dentro de BOX (conservando un punto vecino a cada lado)."""
    inside = [BOX[0] <= x <= BOX[1] and BOX[2] <= y <= BOX[3] for x, y, *_ in line]
    keep = [inside[i] or (i > 0 and inside[i - 1]) or (i < len(line) - 1 and inside[i + 1])
            for i in range(len(line))]
    runs, cur = [], []
    for pt, k in zip(line, keep):
        if k:
            cur.append((round(pt[0], 4), round(pt[1], 4)))
        elif cur:
            runs.append(cur)
            cur = []
    if cur:
        runs.append(cur)
    return [r for r in runs if len(r) >= 2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--download", action="store_true", help="descarga los GeoJSON a data_cache/")
    ap.add_argument("--states-shp", help="shapefile de estados propio (reemplaza a Natural Earth)")
    ap.add_argument("--name-field", default="ENTIDAD", help="campo del .dbf con el nombre del estado")
    a = ap.parse_args()
    CACHE.mkdir(exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    for name, fn in LAYERS.items():
        if name == "states" and a.states_shp:
            shp = Path(a.states_shp)
            attrs = read_dbf(shp.with_suffix(".dbf"))
            geoms = read_shp_polygons(shp)
            assert len(attrs) == len(geoms), "dbf y shp tienen distinto número de registros"
            data = {"features": [
                {"properties": {"name": DISPLAY_NAMES.get(r[a.name_field], r[a.name_field].title()),
                                "iso_a2": "MX"},
                 "geometry": {"type": "Polygon", "coordinates": g}}
                for r, g in zip(attrs, geoms)]}
        else:
            src = CACHE / f"{fn}.geojson"
            if a.download or not src.exists():
                print("Descargando", fn)
                urllib.request.urlretrieve(BASE_URL + fn + ".geojson", src)
            data = json.loads(src.read_text(encoding="utf-8"))
        n = 0
        labels = []
        with open(OUT / f"{name}.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["id", "lon", "lat"])
            for feat in data["features"]:
                p = feat["properties"]
                if name == "states" and p.get("iso_a2") != "MX":
                    continue
                pts = []
                for line in lines_from_geometry(feat["geometry"]):
                    pts += [(x, y) for x, y, *_ in line
                            if DOM[0] <= x <= DOM[1] and DOM[2] <= y <= DOM[3]]
                    for run in clip_runs(line):
                        n += 1
                        for x, y in run:
                            w.writerow([n, x, y])
                if name == "states" and len(pts) >= 15 and p["name"] not in SKIP_LABELS:
                    xs, ys = [q[0] for q in pts], [q[1] for q in pts]
                    lon_c, lat_c = LABEL_POS.get(p["name"], (round((min(xs) + max(xs)) / 2, 3), round((min(ys) + max(ys)) / 2, 3)))
                    labels.append((RENAME.get(p["name"], p["name"]), lon_c, lat_c))
        if labels:
            with open(OUT / "state_labels.csv", "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(["name", "lon", "lat"])
                w.writerows(labels)
        print(f"{name}: {n} polilíneas -> data/basemap/{name}.csv")


if __name__ == "__main__":
    main()
