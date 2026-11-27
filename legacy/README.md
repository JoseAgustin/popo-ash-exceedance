# Script original (`Analisis.R`)

Se conserva sin cambios solo como referencia. Problemas conocidos:

1. `for (i in 1:87600000)` debía ser `1:8760` (o el número real de horas: 8784 en años bisiestos).
2. `Z[1, i, , ]` es un índice incorrecto: con `ncdf4`, `LOAD` llega como `[lon, lat, time, hora]`, no `[time, hora, lat, lon]`.
3. Carga todo `LOAD` en memoria y arma una matriz `L2` de 8760 x 10000 sin necesidad.

La versión corregida está en `R/` y lee por bloques de horas.
