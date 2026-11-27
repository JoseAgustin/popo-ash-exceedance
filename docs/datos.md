# Datos de entrada

Archivos `cXXM-AAAA.nc`, uno por escenario y año, con el depósito de ceniza en suelo del último paso de tiempo de cada pronóstico horario de FALL3D. Los arma `crea_IA.sh` (NCO) a partir de las salidas `p_<HH>h_<CASO>.res.nc` de cada día.

| Escenario | Altura de columna sobre el cráter |
|-----------|-----------------------------------|
| c03M | 3 km |
| c05M | 5 km |
| c10M | 10 km |

El patrón de nombre es `file_pattern: "{scenario}-{year}.nc"` (config). Copias como `c03M-2024.nc.prev_20240424` (versión 7.2 anterior e incompleta de 2024) no coinciden con el patrón y se ignoran. Los años de `years` cuyo archivo no existe se omiten con un aviso.

## Fall3D 7.1 (2018-2023; atributo TITLE "7.2") y Fall3D v8 (2024-2025)

| | 2018-2023 | 2024 y 2025 |
|---|---|---|
| Fall3D | 7.1 (el atributo dice 7.2) | 8 (8.0.1 y 8.2.1) |
| Variable en `c*-AAAA.nc` | `LOAD` | `LOAD` (crea_IA.sh 2.0 lee `tephra_grn_load` y la renombra); se acepta también `tephra_grn_load` |
| Rejilla | 100 × 100 (~0,0404°, ~4,5 km) | 101 × 101 (0,04°) |
| Dominio | lon −100,62 a −96,62; lat 17,02 a 21,02 | el mismo |
| `horas` | escalar igual a 0 (crea_IA.sh 1.0) | un valor por registro: horas desde `AAAA-01-01 00:00` |
| Identificación | atributo global `TITLE = "Fall3d 7.2 results"` (engañoso: la serie es 7.1; el código lo lee tal cual y devuelve 7.2) | atributo `source = "FALL3D model version 8.0.1"` o `Code_version = "8.2.1"` |

Común a ambas versiones: `LOAD(hora, time, lat, lon)` en kg/m2, `time` de tamaño 1, `hora` como dimensión de registro. Cráter: 19,023 N, −98,622 E.

### Versión de Fall3D en 2018-2023
La serie 2018-2023 se calculó con Fall3D 7.1 (confirmado por el autor). El atributo TITLE de los archivos c0XM dice "Fall3d 7.2 results" y el código lo usa tal cual para detectar la versión, por lo que la columna fall3d_version de las tablas muestra 7.2 para 2018-2023: en estas salidas "7.2" significa el valor del atributo, no la versión real.

Los archivos c0XM de abril de 2018 (aprox. 720 horas por escenario) y 8 registros de 2020 se reemplazaron el 2026-11-27 por una reconstrucción desde los .res.nc actuales (título Fall3d 7.1). Una muestra de 174 .res.nc no reemplazados también lleva el título 7.1 y su contenido coincide con el c0XM, coherente con que toda la serie sea 7.1.

Los cinco archivos reemplazados (`c03M-2018`, `c05M-2018`, `c10M-2018`, `c05M-2020`, `c10M-2020`) se documentan en la columna `nota` de `file_inventory.csv` (claves `replaced_files` y `replaced_note` de `config/config.yml`), en el atributo `nota_reemplazo` de `prob_exceedance.nc` y en `docs/correccion_c0XM.md`. Horas reemplazadas: 720 (3 km), 723 (5 km) y 725 (10 km), 1.37-1.38 % de las 52 584 horas de 2018-2023.

## Qué detecta el código en cada archivo
Nada se fija por año; todo se lee del archivo (`python/popo_io.py`, `R/00_utils.R`):

1. **Variable:** la primera de `variable_candidates` (`LOAD`, `tephra_grn_load`) que exista. Si no hay ninguna, error con la lista de variables del archivo.
2. **Dimensiones y rejilla:** número de registros, `lat` y `lon` según las dimensiones reales de la variable (no se asume 100 × 100 ni 8760 registros). El área de celda se calcula con la rejilla de cada archivo (o con la rejilla común, si se interpola).
3. **Versión de Fall3D:** `TITLE` → "7.2" (valor del atributo; la serie 2018-2023 es 7.1); `Code_version` o `source` → "8" (versión completa, p. ej. 8.2.1, en `fall3d_version_full`). Si un archivo v8 no trae esos atributos, se deduce "8" solo si la variable es `tephra_grn_load`; si no, queda "desconocida" con aviso, y se puede fijar con `version_override: ["2024:8"]`.
4. **Tiempo de cada registro:**
   - `horas` escalar o constante (2018-2023, vale 0): el registro `i` (desde 0) es la hora `AAAA-01-01 00:00 + i h`. Esta correspondencia fue **confirmada por el autor de los datos**.
   - `horas` con un valor por registro (v8): se usa ese valor. Si difiere del índice (huecos: crea_IA.sh omite las horas que faltan) se avisa cuántos registros difieren. Si no es entero, creciente y dentro del año, se avisa y se usa el índice.
5. **Valores ausentes:** no finitos, enmascarados o `>= fill_limit` (1e30). Un registro sin ningún valor válido no se cuenta como pronóstico.

## Años completos, parciales y con huecos
El estado de cada año se calcula de las horas con datos (`file_inventory.csv`):

- `completo`: hay datos en las 8760 horas (8784 en bisiesto).
- `parcial`: la última hora con datos es anterior al 31 de diciembre 23:00. La nota indica la fecha de esa hora (p. ej. `parcial (hasta 2025-11-20 13:00)`).
- `con huecos`: llega a fin de año pero faltan horas intermedias (`n_missing_hours`).

Los años parciales no se comparan directamente con años completos (ver `metodologia.md`).

## Inventario verificado (24 archivos, 2018-2025)
| Años | Formato | Variable | Rejilla | Registros | `horas` | Versión (atributo) |
|---|---|---|---|---|---|---|
| 2018, 2019, 2021-2023 | NetCDF clásico 64-bit offset | `LOAD` | 100 × 100 | 8760 | escalar 0 → índice | 7.2 (`TITLE`; real 7.1) |
| 2020 | ídem | `LOAD` | 100 × 100 | 8784 | escalar 0 → índice | 7.2 (`TITLE`; real 7.1) |
| 2024 | NetCDF-4/HDF5 | `LOAD` | 101 × 101 | 8784 | 0..8783 | 8.0.1 (`source`) |
| 2025 | NetCDF-4/HDF5 | `LOAD` | 101 × 101 | 8760 | 0..8759 | 8.2.1 (`Code_version`) |

- Todos los años están completos. En **c05M-2024**, las horas 2024-03-06 02:00 y 03:00 UTC (registros 1562 y 1563, salidas `p_03h_05M` y `p_04h_05M` del 6 de marzo) se calcularon con el WRF del 5 de marzo, porque con el del mismo día Fall3D no resolvió la pluma (`plumeBPT_solve_plume`). El resto de las horas usa el WRF del mismo día.
- **2025 está completo** (8760 h, hasta 2025-12-31 23:00).
- Sin `_FillValue`, sin NaN ni valores de relleno en ningún archivo.
- Los archivos 2018-2023 (TITLE 7.2, real 7.1) tienen valores negativos muy pequeños (ruido numérico, mínimo −0,043 kg/m2); v8 no tiene negativos. No cuentan como excedencia.
- La versión de v8 se toma de los atributos globales del archivo anual, que `ncrcat` copia de la primera salida del año. Si un año mezcla salidas 8.0.1 y 8.2.1 (2024 se completó con 8.0.1), el archivo no lo distingue; en las pruebas de regeneración ambas versiones coincidieron (diferencia < 1e-7 kg/m2 con el mismo WRF).
- Los archivos 7.2 anteriores de 2024 (2290 h, hasta 2024-04-05) fueron sustituidos por los v8.
