# Simulaciones desfasadas en los c0XM (abril de 2018 y 2020): corrección aplicada el 2026-11-27

**Estado: aplicada.** Los cinco archivos afectados (`c03M-2018`, `c05M-2018`, `c10M-2018`, `c05M-2020`, `c10M-2020`) se sustituyeron el 2026-11-27 por versiones corregidas (reemplazo atómico, sha256 antes y después en `reemplazo_MANIFEST.txt`, respaldo de los originales con sha256). `outputs/example/` se regeneró con los datos corregidos.

## Nota sobre la versión de Fall3D

La serie 2018-2023 se calculó con Fall3D 7.1 (confirmado por el autor). El atributo TITLE de los archivos c0XM dice "Fall3d 7.2 results" y el código lo usa tal cual para detectar la versión, por lo que la columna fall3d_version de las tablas muestra 7.2 para 2018-2023: en estas salidas "7.2" significa el valor del atributo, no la versión real.

Los archivos c0XM de abril de 2018 (aprox. 720 horas por escenario) y 8 registros de 2020 se reemplazaron el 2026-11-27 por una reconstrucción desde los .res.nc actuales (título Fall3d 7.1). Una muestra de 174 .res.nc no reemplazados también lleva el título 7.1 y su contenido coincide con el c0XM, coherente con que toda la serie sea 7.1.

El código no cambia: la detección lee `TITLE` y la excepción solo se documenta, en la columna `nota` de `file_inventory.csv` (claves `replaced_files` y `replaced_note` de `config/config.yml`) y en el atributo global `nota_reemplazo` de `prob_exceedance.nc` (transliterado a ASCII, porque el NetCDF clásico no admite acentos en atributos). Los nombres `prob_exc_v72` y la etiqueta `7.2` / `7.2+8` de las tablas se conservan por compatibilidad.

Horas reemplazadas (abril de 2018 más registros de 2020) y fracción de las 52 584 horas de 2018-2023:

| Escenario | Abril 2018 | Registros 2020 | Total reemplazado | Fracción de 2018-2023 |
|---|---|---|---|---|
| 3 km (c03M) | 720 | 0 | 720 | 1.37 % |
| 5 km (c05M) | 720 | 3 | 723 | 1.38 % |
| 10 km (c10M) | 720 | 5 | 725 | 1.38 % |

## Hallazgo

Al comparar el acumulado anual de este repositorio con HazardMaps (mismas simulaciones `p_HHh_CCM.res.nc`, umbral 0.1 kg/m², último paso de cada pronóstico) los conteos por celda coinciden en todos los meses salvo en los siguientes. Las diferencias no vienen del umbral ni del denominador: algunos `.res.nc` se regeneraron **después** de armar los `c0XM` y los registros de los `c0XM` conservaban la versión anterior.

| Archivo c0XM | Registros que difieren | Detalle |
|---|---|---|
| `c03M-2018`, `c05M-2018`, `c10M-2018` | 708, 707 y 706 de 720 en abril de 2018 | Los `.res.nc` de abril de 2018 se regeneraron después de construir el c0XM. Solo coinciden 12-14 registros por caso (las horas 14-24 del 1 de abril y 1-3 más). |
| `c05M-2020` | 3 (2020-08-11 01 h, 2020-10-18 02 h, 2020-11-03 11 h) | `.res.nc` regenerados después de construir el c0XM. |
| `c10M-2020` | 5 (2020-02-03 18 h, 02-06 11 h, 08-03 15 h, 08-07 03 h, 08-15 21 h) | Los registros del c0XM tenían un `LOAD` máximo de ~2e-4 kg/m2; los `.res.nc` actuales dan 0.5-3.4 kg/m2. |

En el resto de 2018-2023 no se encontraron diferencias: se compararon bit a bit los 3 850 `.res.nc` con fecha de modificación posterior a la construcción del `c0XM` (menos 24 h) y 720 simulaciones de control al azar (40 por año y caso); todas las demás coinciden. Dos celdas en 10M 2019 (febrero y junio) difieren solo porque HazardMaps usa `>=` y este repositorio `>`: el valor es exactamente float32(0.1).

**Diferencias entre las dos corridas de abril de 2018.** Los `.res.nc` actuales (título Fall3d 7.1) tienen en abril de 2018 una masa depositada de 0.35-0.48 veces la de los `c0XM` anteriores (mediana por simulación, según el caso); la correlación espacial mediana del campo diario medio es 0.89 en 3 km y 0.77 en 10 km. El autor eligió la reconstrucción desde los `.res.nc` actuales (la misma que usó HazardMaps para el artículo).

## Procedimiento

`scripts/corrige_registros_c0XM.py` copia el original a un archivo nuevo (los `.res.nc` leídos del servidor en solo lectura) y sustituye solo los registros indicados por el `LOAD` del último paso de `time` del `.res.nc` (equivale a `ncks -d time,-1 -v LOAD,lat,lon`), con el registro `i = AAAA-01-01 00:00 + i horas`. No cambia atributos, formato (NetCDF3 de 64 bits), dimensiones ni otras variables, y verifica bit a bit el resultado.

```
python scripts/corrige_registros_c0XM.py --original c03M-2018.nc --reemplazos rep_c03M-2018.npz --salida c03M-2018.corregido.nc
```

## Efecto sobre los resultados (corrida con el mismo código, 18 archivos 2018-2023 y 6 v8)

Con los `c0XM` corregidos, `prob_exc_year` coincide con HazardMaps con diferencia 0 en 2018-2022, salvo las 2 celdas de 10M 2019. Frente a los valores actuales de `outputs/` (que la corrida de referencia reproduce byte a byte):

| Magnitud | 2018 | 2018-2023 | 2018-2025 |
|---|---|---|---|
| P máxima, 3 km (%) | 16.735 → 16.632 | 16.729 → 16.726 | 15.687 → 15.674 |
| P máxima, 5 km (%) | 17.603 → 17.591 | 17.117 → 17.112 | sin cambio |
| P máxima, 10 km (%) | 15.833 → 15.799 | 14.449 → 14.503 | 13.785 → 13.773 |
| Cambio máximo del área por nivel de probabilidad | -7.7 % (3 km, 5 %) | -1.3 % | -1.4 % |
| Horas con excedencia (3 km / 10 km) | -3 / sin cambio | -3 / +5 | -3 / +5 |

La posición de la probabilidad máxima no cambia, salvo en 10M 2018, donde la longitud pasa de -97.711 a -97.630 (un nodo de la rejilla). Los periodos y las conclusiones no cambian.

## Pasos realizados
1. Decisión del usuario (opción A): sustituir los cinco archivos.
2. Respaldo de los originales (sha256) y reemplazo atómico por archivo (copia temporal, verificación, `mv`).
3. Reejecución de `python/exceedance.py` y `R/run_all.R` (NetCDF idénticos, tablas dentro de 1e-6), test sintético y regeneración de `outputs/example/`.
4. Los cambios en tablas y figuras de `outputs/example/` respecto al repositorio publicado se describen en el historial de git.

## Pendientes (no realizados; decide el usuario)
Las salidas aún no reflejan el texto actual de la nota (`outputs/example/` y `outputs/` se generaron con la nota anterior, que hablaba de una excepción de 7.1 dentro de una serie 7.2). Regenerarlas exige volver a correr `python/exceedance.py` y `R/run_all.R`, y copiar de nuevo los archivos de `outputs/example/` que contienen la nota o el valor `7.2`: `file_inventory.csv`, `max_probability.csv`, `summary_scenario_period.csv`, `area_by_probability_level.csv` y `period_definitions.csv`. Los PNG llevan el texto "Fall3D 7.2" en sus figuras.
