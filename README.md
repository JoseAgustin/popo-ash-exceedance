# popo-ash-exceedance

Análisis de excedencia del depósito de ceniza (`LOAD`, kg/m2) en pronósticos horarios FALL3D del volcán Popocatépetl, para columnas de 3, 5 y 10 km (c03M, c05M, c10M), 2018-2025.

| Años | Fall3D | Rejilla | Estado |
|---|---|---|---|
| 2018-2023 | 7.2 | 100 × 100 | completos |
| 2024 | 8 (archivo: 8.0.1) | 101 × 101, mismo dominio | completo (8784 h); en c05M, 2024-03-06 02:00 y 03:00 UTC se calcularon con el WRF del día previo |
| 2025 | 8 (archivo: 8.2.1) | 101 × 101 | completo (8760 h, hasta 2025-12-31 23:00) |

El código detecta en cada archivo la variable (`LOAD` o `tephra_grn_load`), la rejilla, la versión de Fall3D y el tiempo de cada registro (`horas` en 2024-2025), y marca los años parciales o con huecos a partir de sus horas reales.

Periodos que se reportan:
- **cada año**, en la rejilla común 100 × 100 (2024 y 2025 interpolados de 101 × 101);
- **2018-2023**: solo Fall3D 7.2, rejilla nativa 100 × 100 (referencia homogénea);
- **2024-2025**: solo Fall3D v8, rejilla nativa 101 × 101, **sin interpolar**;
- **2018-2025**: todos los años en la rejilla común, con advertencia de mezcla de versiones.

**Una diferencia entre 2018-2023 y 2024-2025 puede reflejar el cambio de versión del modelo (7.2 → 8) y no solo la meteorología.** Detalles en `docs/datos.md` y `docs/metodologia.md`.

Umbral por defecto: `LOAD > 0.1 kg/m2` (parámetro `threshold` en `config/config.yml`).

## Estructura
```
config/config.yml   parámetros (ruta de datos, umbral, escenarios, años, rejilla)
R/                  implementación principal (run_all.R es el punto de entrada)
python/             implementación equivalente (exceedance.py; lectura y rejilla en popo_io.py)
LICENSE, NOTICE     MIT para el código; condiciones de los datos de terceros (CONABIO, Natural Earth)
tests/              test sintético de equivalencia Python/R (no necesita los datos reales)
data/basemap/       mapa base recortado en CSV: estados (CONABIO), costa (Natural Earth), ciudades
scripts/            make_basemap.py (regenera data/basemap)
docs/               metodología y descripción de los datos
legacy/             script original y sus errores
outputs/            resultados (ignorados por git, salvo outputs/example/)
```

## Uso
Edita `data_dir` en `config/config.yml` (los `.nc` no se incluyen en el repositorio).

**R** (>= 4.0; paquetes `ncdf4`, `yaml`, `data.table`, `ggplot2`):
```
Rscript -e 'install.packages(c("ncdf4","yaml","data.table","ggplot2"))'
Rscript R/run_all.R config/config.yml [carpeta_de_salida]
```
Si no se indica carpeta, se usa `output_dir` de la configuración. En macOS, `install.packages` usa binarios de CRAN y no requiere instalar libnetcdf.

**Python** (>= 3.9):
```
pip install -r python/requirements.txt
python python/exceedance.py --config config/config.yml
# prueba rápida
python python/exceedance.py --scenarios c03M --years 2018 --output-dir outputs_test
```
Con R o Python, 18 archivos de un año (2018-2023) se procesan en ~1 minuto y con poca memoria (bloques de 1000 horas).

Python lee NetCDF clásico con scipy; si los archivos son NetCDF-4/HDF5 hace falta el paquete `netCDF4` (el script lo indica). R (`ncdf4`) lee ambos formatos.

### Opciones de rejilla (`config/config.yml`)
- `grid_policy: regrid` y `regrid_n: 100` (**decisión del proyecto**, fijada en `config.yml`): cada año y el combinado de todos en 100 × 100; 2018-2023 sin interpolar, 2024-2025 interpolados. Los combinados de una sola versión (2018-2023 y 2024-2025) se calculan siempre en su rejilla nativa. `regrid_n: 101` llevaría 2018-2023 a la rejilla v8.
- `grid_policy: native`: cada año en su rejilla; no hay combinado de todos los años y se escribe un NetCDF por rejilla.
- También por línea de comandos en Python: `--grid-policy native`.

## Entorno de Python
Probado con Python 3.11.18 en un entorno virtual creado fuera del repositorio (numpy 2.4.6, scipy 1.17.1, pandas 3.1.1, matplotlib 3.11.2, netCDF4 1.7.5.1, PyYAML 6.0.4; todos instalados desde binarios). También funciona con Python 3.9 (numpy 2.0, pandas 2.3).

Crear el entorno (una vez):
```
python3.11 -m venv /ruta/a/venv_popo
source /ruta/a/venv_popo/bin/activate
pip install --only-binary=:all: -r python/requirements.txt
```
Usarlo:
```
source /ruta/a/venv_popo/bin/activate
python python/exceedance.py --config config/config.yml
deactivate
```
- Necesarios: `numpy`, `scipy`, `pandas`, `matplotlib`.
- Recomendados: `netCDF4` (lee NetCDF-4/HDF5 y clásico; si falta, `scipy` lee solo el formato clásico) y `PyYAML` (si falta, se usa un lector mínimo de `config.yml`). Si `netCDF4` está instalado, el script lo usa en lugar de `scipy`.
- `--only-binary=:all:` evita compilar paquetes; si para alguna versión de Python no hay binarios, pip falla en vez de compilar y conviene usar otra versión (p. ej. 3.11).
- Un venv queda ligado al intérprete con que se creó: si el gestor de paquetes elimina esa versión de Python (p. ej. al actualizar Homebrew), hay que recrear el venv.
- Comprobación sin datos reales: `tests/run_synthetic_tests.sh /tmp/popo_synth /ruta/a/venv_popo/bin/python` (requiere también R con `ncdf4`).

## Test sintético
Crea archivos pequeños tipo 7.2 y v8 fuera del repositorio, corre Python y R con `regrid` y `native` y compara sus salidas:
```
tests/run_synthetic_tests.sh /tmp/popo_synth [python3]
```

## Salidas (`outputs/`)
- `tables/summary_scenario_period.csv`: por escenario, año y combinado (horas con excedencia, área en km2, percentiles).
- `tables/summary_monthly.csv`, `summary_seasonal.csv`, `summary_diurnal.csv`.
- `tables/max_probability.csv`, `area_by_probability_level.csv`.
- `tables/exceedance_hours.csv`, `hourly_metrics.csv.gz` (con la hora real de cada registro).
- `tables/file_inventory.csv`: por archivo, variable, versión de Fall3D, rejilla nativa y de trabajo, horas con datos, primera y última hora, estado (`completo`, `parcial`, `con huecos`) y nota.
- `tables/period_definitions.csv`: años de cada periodo, versión, rejilla y si se combinó.
- Todas las tablas por periodo llevan `fall3d_version`, `grid`, `grid_native` y `nota` (año parcial, mezcla de versiones).
- Periodos: cada año, 2018-2023 (solo 7.2), 2024-2025 (solo v8, rejilla nativa 101 × 101) y 2018-2025 (7.2+8).
- `netcdf/prob_exceedance.nc` (rejilla común 100 × 100): `prob_exc_year(scenario, year, lat, lon)` de 2018-2025, `prob_exc_all` (2018-2025), `prob_exc_v72` (2018-2023), conteos (`count_exc_*`) y horas válidas (`n_valid_*`).
- `netcdf/prob_exceedance_101x101.nc` (rejilla nativa v8, sin interpolar): `prob_exc_year` de 2024 y 2025 y `prob_exc_v8` (2024-2025).
- Atributos globales: `grid`, `grid_space`, `fall3d_version_by_year`, `partial_years`, `years_with_gaps` (p. ej. `2024:2h`).
- `figures/*.png`: mapas combinados 2018-2023 (`maps_prob_combined.png`), 2024-2025 (`maps_prob_combined_v8.png`) y 2018-2025 (`maps_prob_combined_all.png`); por año (con versión, interpolación y huecos), diferencias entre escenarios (2018-2023), barras por año, área y estacionalidad.
- `outputs/example/` (versionado): copia pequeña de mapas combinados, diferencias y tablas resumen.

## Estado de verificación
- Datos reales 2018-2025 (24 archivos): R y Python dan NetCDF idénticos (probabilidades, conteos y horas válidas, diferencia 0, en ambas rejillas) y tablas con diferencia relativa < 1e-7 (máximo de LOAD en float32). Los resultados por año 2018-2023 no cambian respecto a versiones anteriores del código.
- Test sintético: 155 controles correctos con Python y R.
- Versiones usadas: R 4.5.3, ncdf4 1.24, yaml 2.3.12, data.table 1.18.2, ggplot2 4.0.2; Python 3.11.18 (venv, con netCDF4) y Python 3.9.6 (numpy 2.0.2, scipy 1.13.1, pandas 2.3.3, matplotlib 3.9.4, sin netCDF4). Ambos Python dan los mismos resultados.

Tiempo de los registros: en 2018-2023 el registro `i` corresponde a `AAAA-01-01 00:00 + i horas` (confirmado por el autor de los datos); en 2024-2025 se usa la variable `horas` de cada registro (ver `docs/datos.md`).

Ver `docs/metodologia.md` para definiciones y limitaciones.

## Licencia
- **Código**: MIT (`LICENSE`).
- **Datos del mapa base** (`data/basemap/`): **no** están bajo MIT y conservan sus propias condiciones; ver `NOTICE`. Los estados (`states.csv`, `state_labels.csv`) son de CONABIO (2005), *División Política Estatal* 1:1 000 000, con información del INEGI (2000): **uso no lucrativo y se debe citar la fuente del mapa y a CONABIO**. La costa es de Natural Earth (dominio público).

## Mapa base y atribución
- **Estados:** CONABIO (2005), *División Política Estatal 1:1 000 000* (`destdv1gw`), elaborada por CONABIO con información cartográfica del INEGI (Conjunto de datos vectoriales topográficos y toponímicos 1:1 000 000, 2000). Coordenadas geográficas WGS84. Según los metadatos de CONABIO, **no se permite utilizar estos datos con fines lucrativos y se debe citar la fuente del mapa y a CONABIO**. Fuente: <http://www.conabio.gob.mx/informacion/gis/>.
- **Costa:** [Natural Earth](https://www.naturalearthdata.com/) 1:10m (dominio público), del espejo [nvkelso/natural-earth-vector](https://github.com/nvkelso/natural-earth-vector).

El repositorio incluye solo polilíneas recortadas al dominio (`data/basemap/*.csv`, ~0.4 MB). La licencia MIT cubre el código; los datos de estados (`states.csv`, `state_labels.csv`) conservan las condiciones de CONABIO indicadas arriba (cita completa en `NOTICE`). Las posiciones de ciudades y de etiquetas de estados se definieron manualmente.

Para regenerar el mapa base (solo biblioteca estándar de Python):
```
# costa de Natural Earth (descarga ~50 MB a data_cache/, ignorado por git) y estados desde un shapefile propio
python scripts/make_basemap.py --states-shp RUTA/destdv1gw.shp --name-field ENTIDAD
# alternativa: estados de Natural Earth
python scripts/make_basemap.py --download
```
Los scripts de R y Python leen `data/basemap/` desde la raíz del repositorio. Ejecuta `Rscript R/run_all.R` desde la raíz.

## Autoría y contacto
Dr. José Agustín García Reynoso, Instituto de Ciencias de la Atmósfera y Cambio Climático (ICAyCC), Universidad Nacional Autónoma de México (UNAM). Correo: agustin@atmosfera.unam.mx

Cita sugerida (sin DOI):
> García Reynoso, J. A. (2026). *popo-ash-exceedance: análisis de excedencia de depósito de ceniza en pronósticos FALL3D del volcán Popocatépetl* (código R y Python). ICAyCC, UNAM.

## Ruta de los datos
`data_dir` en `config/config.yml` es un valor de ejemplo (`/ruta/a/los/netcdf`). Prioridad: `--data-dir` (Python) > variable de entorno `POPO_DATA_DIR` > `config.yml`. Para uso local puedes crear `config/config.local.yml` (ignorado por git):
```
POPO_DATA_DIR=/mi/carpeta/con/nc Rscript R/run_all.R config/config.yml
POPO_DATA_DIR=/mi/carpeta/con/nc python python/exceedance.py
```
