"""Lectura robusta de los NetCDF cXXM-AAAA.nc (Fall3D 7.1/7.2 y v8) y utilidades de rejilla.

Autor: Dr. José Agustín García Reynoso (ICAyCC, UNAM) <agustin@atmosfera.unam.mx>
Licencia: MIT (ver LICENSE).

Qué se detecta en cada archivo (nada se asume por el año):
  * variable de depósito: la primera de `variable_candidates` presente (LOAD, tephra_grn_load);
  * dimensiones (registros, lat, lon) a partir de las dimensiones reales de esa variable;
  * versión de Fall3D: atributo TITLE ("Fall3d 7.2 results"; en 2018-2023 el atributo es engañoso: la serie es Fall3D 7.1) o Code_version / source ("... 8.0.1");
  * tiempo de cada registro: variable `horas` si trae un valor por registro; si es escalar o
    constante (caso 2018-2023, vale 0) se usa el índice del registro (hora del año).
"""
import re
import warnings
from pathlib import Path

import numpy as np

EARTH_RADIUS_KM = 6371.0088


# --------------------------------------------------------------------------
# Versión del modelo y tiempo
# --------------------------------------------------------------------------
def detect_version(attrs, varname, override=None):
    """Devuelve (etiqueta, versión_completa, origen). Etiqueta: '7.2' (valor del atributo TITLE; 2018-2023 es en realidad 7.1), '8' o 'desconocida'."""
    if override:
        full = str(override)
        return (full if not full.startswith("8") else "8"), full, "config:version_override"
    attrs = {str(k): (v.decode("utf-8", "replace") if isinstance(v, bytes) else v) for k, v in attrs.items()}
    title = str(attrs.get("TITLE", ""))
    m = re.search(r"fall3d\s*v?(\d+\.\d+(?:\.\d+)?)", title, re.I)
    if m:
        full = m.group(1)
        return _label(full), full, "TITLE"
    for key in ("Code_version", "code_version", "source"):
        if key in attrs:
            m = re.search(r"(\d+\.\d+(?:\.\d+)?)", str(attrs[key]))
            if m:
                full = m.group(1)
                return _label(full), full, key
    if varname == "tephra_grn_load":
        return "8", "8", "variable"
    return "desconocida", "desconocida", "ninguno"


def _label(full):
    """'7.2.1' -> '7.2'; '8.0.1' -> '8' (en v8 solo se distingue la versión mayor)."""
    parts = full.split(".")
    return parts[0] if int(parts[0]) >= 8 else f"{parts[0]}.{parts[1]}"


def expected_hours(year):
    return 8784 if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)) else 8760


def resolve_hours(horas, nrec, year):
    """Hora del año (desde 0) de cada registro. Devuelve (idx[int64], fuente, n_discrepancias, avisos)."""
    warns = []
    idx = np.arange(nrec, dtype=np.int64)
    if horas is None:
        return idx, "index", 0, warns
    h = np.asarray(horas).ravel()
    if h.size <= 1 or (h.size == nrec and np.all(h == h[0])):
        return idx, "index", 0, warns          # escalar o constante (p. ej. 0): vale el índice
    if h.size != nrec:
        warns.append(f"'horas' tiene {h.size} valores y hay {nrec} registros; se usa el índice")
        return idx, "index", 0, warns
    hi = np.rint(h).astype(np.int64)
    if (np.abs(h - hi) > 1e-6).any() or np.any(np.diff(hi) <= 0) or hi[0] < 0 or hi[-1] >= expected_hours(year):
        warns.append("'horas' no es entera, creciente y dentro del año; se usa el índice")
        return idx, "index", 0, warns
    n_mis = int((hi != idx).sum())
    if n_mis:
        warns.append(f"'horas' difiere del índice en {n_mis} registros (huecos en la serie); se usa 'horas'")
    return hi, "horas", n_mis, warns


# --------------------------------------------------------------------------
# Lectura
# --------------------------------------------------------------------------
class NCFile:
    """Acceso por bloques a la variable de depósito (registro, [dims de tamaño 1], lat, lon)."""

    def __init__(self, path, var_candidates=("LOAD", "tephra_grn_load")):
        self.path = str(path)
        with open(self.path, "rb") as fh:
            is_hdf5 = fh.read(8).startswith(b"\x89HDF")
        try:
            import netCDF4  # type: ignore
            self._kind = "netCDF4"
            self._ds = netCDF4.Dataset(self.path)
            self._vars = self._ds.variables
            self.attrs = {k: self._ds.getncattr(k) for k in self._ds.ncattrs()}
        except ImportError:
            if is_hdf5:
                raise RuntimeError(f"{self.path} es NetCDF-4/HDF5 y falta el paquete netCDF4 "
                                   "(pip install netCDF4); alternativa: usar la versión en R (ncdf4).")
            from scipy.io import netcdf_file
            self._kind = "scipy"
            self._ds = netcdf_file(self.path, mmap=True)
            self._vars = self._ds.variables
            self.attrs = dict(self._ds._attributes)
        self.varname = next((v for v in var_candidates if v in self._vars), None)
        if self.varname is None:
            raise KeyError(f"{Path(self.path).name}: no hay ninguna de las variables {list(var_candidates)}; "
                           f"hay {list(self._vars)}")
        self._var = self._vars[self.varname]
        dims = tuple(self._var.dimensions)
        dl = [str(d).lower() for d in dims]
        try:
            self._ilat = next(i for i, d in enumerate(dl) if d.startswith("lat"))
            self._ilon = next(i for i, d in enumerate(dl) if d.startswith("lon"))
        except StopIteration:
            raise KeyError(f"{Path(self.path).name}: {self.varname}{dims} no tiene dimensiones lat/lon")
        if 0 in (self._ilat, self._ilon):
            raise KeyError(f"{Path(self.path).name}: la primera dimensión de {self.varname} debe ser el registro")
        shape = tuple(self._var.shape)
        self.nrec, self.nlat, self.nlon = int(shape[0]), int(shape[self._ilat]), int(shape[self._ilon])
        self.dims = dims
        self.lat = self._get(dims[self._ilat]).astype("f8")
        self.lon = self._get(dims[self._ilon]).astype("f8")
        self.horas = self._get("horas") if "horas" in self._vars else None

    def _get(self, name):
        v = self._vars[name]
        if self._kind == "scipy":
            return np.array(v.data)
        return np.ma.filled(np.asarray(v[...]), np.nan)

    def read(self, h0, h1):
        """float32 (n, lat, lon); relleno enmascarado -> NaN. Dimensiones extra (time) -> último valor."""
        sl = []
        for i in range(len(self.dims)):
            sl.append(slice(h0, h1) if i == 0 else slice(None) if i in (self._ilat, self._ilon) else -1)
        blk = self._var[tuple(sl)]
        blk = np.ma.filled(blk.astype("f4"), np.nan) if self._kind == "netCDF4" else np.array(blk, dtype="f4")
        if self._ilat > self._ilon:
            blk = np.swapaxes(blk, 1, 2)
        return blk

    def close(self):
        self._var = None
        self._vars = None
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                self._ds.close()
            except Exception:
                pass


# --------------------------------------------------------------------------
# Rejilla
# --------------------------------------------------------------------------
def cell_area_km2(lat, lon):
    """Área de cada celda (km2) en malla regular lat/lon: R^2 dlon dlat cos(lat). Matriz [lat, lon]."""
    dlat = np.deg2rad(np.abs(np.mean(np.diff(lat))))
    dlon = np.deg2rad(np.abs(np.mean(np.diff(lon))))
    a = EARTH_RADIUS_KM ** 2 * dlat * dlon * np.cos(np.deg2rad(lat))
    return np.repeat(a[:, None], len(lon), axis=1)


def target_axes(cfg):
    n = int(cfg["regrid_n"])
    return (np.linspace(cfg["domain_lat"][0], cfg["domain_lat"][1], n),
            np.linspace(cfg["domain_lon"][0], cfg["domain_lon"][1], n))


def is_same_grid(lat, lon, tlat, tlon, tol=1e-4):
    return (len(lat) == len(tlat) and len(lon) == len(tlon)
            and np.max(np.abs(lat - tlat)) < tol and np.max(np.abs(lon - tlon)) < tol)


def axis_weights(src, tgt, tol=1e-4):
    """Índices (0-based) del nodo inferior, peso del superior y marca 'fuera del dominio'."""
    src = np.asarray(src, "f8")
    tgt = np.asarray(tgt, "f8")
    if np.any(np.diff(src) <= 0):
        raise ValueError("las coordenadas de la rejilla deben ser estrictamente crecientes")
    i0 = np.clip(np.searchsorted(src, tgt, side="right") - 1, 0, len(src) - 2)
    w = np.clip((tgt - src[i0]) / (src[i0 + 1] - src[i0]), 0.0, 1.0)
    outside = (tgt < src[0] - tol) | (tgt > src[-1] + tol)
    return i0, w, outside


class Regridder:
    """Interpolación bilineal (en doble precisión) de (n, lat, lon) a otra rejilla regular.

    Primero a lo largo de lat y luego de lon; un NaN en cualquiera de las 4 esquinas da NaN.
    El resultado se redondea a float32 (como los datos) antes de compararlo con el umbral.
    """

    def __init__(self, src_lat, src_lon, tgt_lat, tgt_lon):
        self.iy, self.wy, oy = axis_weights(src_lat, tgt_lat)
        self.ix, self.wx, ox = axis_weights(src_lon, tgt_lon)
        self.outside = oy[:, None] | ox[None, :]

    def __call__(self, a):
        a = np.asarray(a, dtype="f8")
        wy = self.wy[None, :, None]
        b = (1.0 - wy) * a[:, self.iy, :] + wy * a[:, self.iy + 1, :]
        wx = self.wx[None, None, :]
        c = (1.0 - wx) * b[:, :, self.ix] + wx * b[:, :, self.ix + 1]
        c[:, self.outside] = np.nan
        return c.astype("f4")
