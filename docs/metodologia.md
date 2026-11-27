# Metodología

## Excedencia
Para cada celda y cada conjunto de horas: `P = N(LOAD > u) / N(valores válidos)`, con umbral `u` (por defecto 0.1 kg/m2, parámetro `threshold`) y comparación estricta. NA y relleno se excluyen del denominador.

Como cada hora es un pronóstico independiente de una erupción hipotética (escenario de columna fija), `P` es la fracción de esos pronósticos en que el depósito supera `u`. **No es la probabilidad de que ocurra una erupción.** Es una probabilidad condicional a que ocurra el evento simulado, bajo la meteorología de cada hora.

## Rejillas distintas: `grid_policy`
Fall3D 7.1 (2018-2023; el atributo TITLE de los archivos dice 7.2) usa 100 × 100 nodos y Fall3D v8 (2024-2025) 101 × 101 sobre el mismo dominio. Las celdas no se mezclan directamente.

- **`regrid`** (decisión del proyecto, con `regrid_n: 100`): `LOAD` de cada archivo cuya rejilla no coincide con la rejilla común se interpola de forma bilineal a esa rejilla, **antes** de compararlo con el umbral. La rejilla común tiene `regrid_n` × `regrid_n` nodos sobre `domain_lon`/`domain_lat`: `100` es la de 2018-2023 (esos años quedan sin interpolar) y `101` la de v8. Se usa para cada año y para el combinado de todos los años. Los combinados de una sola versión se calculan en su rejilla nativa (ver "Periodos").
- **`native`**: cada año se analiza en su rejilla. Solo se combinan años con la misma rejilla; un periodo con rejillas distintas se omite (con aviso y `status = omitido` en `period_definitions.csv`). Hay un NetCDF por rejilla (`prob_exceedance_100x100.nc`, `prob_exceedance_101x101.nc`).

Las tablas registran la rejilla de cada fila: `grid` (rejilla de trabajo) y `grid_native` (rejilla del archivo); `file_inventory.csv` añade `regridded`.

### Por qué interpolar LOAD y no la probabilidad o la excedencia
La excedencia es una función no lineal (escalón) de `LOAD`. Interpolar el campo continuo y después aplicar el umbral conserva el significado de "depósito mayor que 0,1 kg/m2 en ese nodo" en cada pronóstico; interpolar el indicador 0/1 o la probabilidad ya calculada daría valores intermedios que no corresponden a ningún pronóstico y no permitiría sumar conteos entre años.

### Sesgo posible
- La interpolación bilineal suaviza: no crea máximos nuevos y puede rebajar picos estrechos (de una o dos celdas). Cerca del umbral, una celda interpolada puede quedar por debajo de 0,1 kg/m2 aunque un nodo vecino lo supere, o al revés. El efecto se espera pequeño porque las dos rejillas casi coinciden (0,0404° frente a 0,04°; desplazamiento máximo ~0,02°, la mitad de una celda), pero concentrado en el borde de la pluma y en c03M, cuyo depósito es más estrecho.
- `max_load_kg_m2` de los años interpolados es el máximo del campo interpolado (puede ser algo menor que el del archivo).
- Medido con los datos reales 2024-2025 (interpolado a 100 × 100 frente a nativo 101 × 101):

| Escenario | Horas con excedencia | Área media (todas las horas) | P máxima | LOAD máximo |
|---|---|---|---|---|
| c03M | 16 206 frente a 16 669 (−2,8 %) | 274 frente a 299 km2 (−8,5 %) | 0,137 frente a 0,138 | 0,73 frente a 1,12 kg/m2 |
| c05M | igual (17 542) | 1905 frente a 1870 km2 (+1,9 %) | 0,224 frente a 0,216 | 2,81 frente a 3,88 kg/m2 |
| c10M | igual (17 544) | 5054 frente a 4895 km2 (+3,2 %) | 0,179 frente a 0,173 | 11,87 frente a 13,50 kg/m2 |

  El suavizado recorta los picos (LOAD máximo 12-35 % menor) y, en c03M, las plumas estrechas que apenas superan el umbral. En c05M y c10M el área crece un poco, en parte porque las celdas de 100 × 100 son ~2 % más grandes. Por eso el combinado 2024-2025 se calcula sin interpolar y los años v8 en la rejilla común deben leerse con este sesgo en mente.

## Periodos
- Cada año por separado (`period = AAAA`).
- **Combinados homogéneos** (`kind = homogeneo`): uno por versión de Fall3D con al menos dos años, calculado en la **rejilla nativa** de esa versión, sin interpolar: 2018-2023 (7.1, etiqueta `7.2` por el atributo TITLE; 100 × 100; ver `docs/correccion_c0XM.md`) y 2024-2025 (v8, 101 × 101). Son las referencias para comparar dentro de una versión. En el NetCDF: `prob_exc_v72` (en `prob_exceedance.nc`) y `prob_exc_v8` (en `prob_exceedance_101x101.nc`).
- **Combinado de todos los años** (`kind = todos`): 2018-2025 (etiqueta `7.2+8`; real 7.1 + 8) en la rejilla común. Su nota advierte que una diferencia puede deberse al cambio de versión del modelo.
- Si todos los años son de la misma versión solo hay un combinado (`todos`). `period_definitions.csv` lista los periodos con años, versión, rejilla y estado.
- 2018-2023 y 2024-2025 tienen rejillas distintas, así que sus mapas no se restan celda a celda; para compararlos usa las tablas (área, horas, probabilidad máxima) o los años v8 en la rejilla común.

Los combinados suman conteos de excedencia y horas válidas de sus años, de modo que cada año pesa según sus horas. Un año parcial pesa menos y su estación incompleta sesga el combinado; la nota lo indica (`incluye año parcial 2025 (hasta ...)`).

## Cambio de versión del modelo (7.1 → 8; etiqueta 7.2 → 8)
La física, la parametrización y la configuración de Fall3D 7.1 y 8 no son idénticas, y la rejilla cambia. Una diferencia entre 2018-2023 y 2024-2025 puede reflejar el cambio de versión del modelo y no solo la variabilidad meteorológica entre años. Por eso se reporta el combinado homogéneo 2018-2023 por separado, las tablas llevan la columna `fall3d_version` y las figuras indican la versión de cada año. Para atribuir diferencias al año haría falta correr un mismo periodo con ambas versiones.

## Años parciales
El estado (`completo`, `parcial`, `con huecos`) se calcula de las horas con datos de cada archivo, no se fija en la configuración (ver `datos.md`). Un año parcial no cubre todas las estaciones, así que su fracción de horas con excedencia y su mapa de probabilidad **no se comparan directamente** con años completos; para compararlo, usa los mismos meses (`summary_monthly.csv`). Las tablas (`nota`), el NetCDF (`partial_years`, `years_with_gaps`) y las figuras (título, asterisco y achurado) marcan los años parciales con la fecha de su última hora y los años con huecos con las horas faltantes. En los datos actuales no hay años parciales ni con huecos.

## Métricas horarias
- `ncells`: celdas con LOAD > u en esa hora. `area_km2`: suma de áreas de esas celdas, con la rejilla de trabajo.
- Área de celda: `R^2 * dlat * dlon * cos(lat)` (R = 6371.0088 km, radianes).
- `time` (en `exceedance_hours.csv` y `hourly_metrics.csv.gz`): hora real de cada registro (`horas` en v8; índice en 2018-2023).
- "Hora con excedencia": hora con al menos una celda sobre el umbral. Esta métrica casi siempre vale ~100% (hay depósito cerca del cráter en casi todos los pronósticos); la información útil está en el mapa de probabilidad y en el área.
- Percentiles de área (`p50`, `p90`, `p99`) sobre horas con excedencia.
- `area_by_probability_level.csv`: área de celdas con `P >= nivel` (niveles en `prob_levels`).

## Comparación con el umbral (float32)
`LOAD` se guarda en float32. El umbral se redondea a float32 antes de comparar (0.1 pasa a 0.100000001490116) para que un valor almacenado igual a 0.1f no cuente como excedencia. La interpolación se calcula en doble precisión y el resultado se redondea a float32 antes del umbral, igual en R y en Python.

## Comparación R vs Python
- Datos reales 2018-2025 (24 archivos, rejilla común y nativa v8): NetCDF idénticos (diferencia 0 en probabilidades, conteos y horas válidas) y tablas con diferencia relativa < 1e-7 (máximo de LOAD, conversión float32/double).
- Test sintético (`tests/run_synthetic_tests.sh`): archivos con TITLE 7.2 (100 × 100, `horas` = 0) y v8 (101 × 101, `tephra_grn_load` y `LOAD`, `horas` con huecos, un registro de relleno, años completos y parcial), con `regrid` y `native`. Python y R dan los mismos conteos por celda, NetCDF y tablas (diferencia relativa < 1e-7); 155 controles. La interpolación reproduce sin error un campo lineal.
- Las configuraciones YAML con `1.0e30` se leen como texto; usar `1.0e+30`.

## Limitaciones
- El dominio mide 4° x 4° (~190 000 km2). Con columna de 10 km el depósito llega a los bordes, así que áreas y probabilidades de c10M están truncadas por el dominio.
- Los mapas no tienen proyección: se dibujan en lon/lat con relación de aspecto 1/cos(19°). Estados: CONABIO/INEGI 1:1 000 000 (uso no lucrativo, citar fuente); costa: Natural Earth 1:10m (dominio público). La costa de Natural Earth y el límite del shapefile difieren hasta ~3 km en el dominio (menos de una celda).
- Estaciones: lluvias = mayo-octubre, secas = noviembre-abril.
