#!/bin/bash
# Test de equivalencia Python/R con datos sintéticos (Fall3D 7.2 100x100 y v8 101x101), sin datos reales.
# Autor: Dr. José Agustín García Reynoso (ICAyCC, UNAM) <agustin@atmosfera.unam.mx>. Licencia: MIT.
# Uso (desde la raíz del repositorio):
#   tests/run_synthetic_tests.sh [carpeta_temporal] [python]
# La carpeta temporal (por omisión /tmp/popo_synth) queda fuera del repositorio.
set -u
TMP=${1:-/tmp/popo_synth}
PY=${2:-python3}
case "$(cd "$(dirname "$TMP")" 2>/dev/null && pwd)/$(basename "$TMP")" in
  "$(pwd)"/*) echo "ERROR: la carpeta temporal debe estar fuera del repositorio"; exit 2 ;;
esac
rm -rf "$TMP"
Rscript tests/make_synthetic_data.R "$TMP" || exit 2
status=0
for c in dataset1_regrid dataset1_native dataset2; do
  echo "=== Caso $c"
  "$PY" python/exceedance.py --config "$TMP/$c.yml" --output-dir "$TMP/out_${c}_Py" > "$TMP/$c.py.log" 2>&1 || { echo "Python falló (ver $TMP/$c.py.log)"; status=1; continue; }
  Rscript R/run_all.R "$TMP/$c.yml" "$TMP/out_${c}_R" > "$TMP/$c.R.log" 2>&1 || { echo "R falló (ver $TMP/$c.R.log)"; status=1; continue; }
  "$PY" tests/compare_outputs.py "$TMP/out_${c}_Py" "$TMP/out_${c}_R" --expect "$c" || status=1
done
echo "=== Caso nc4 (NetCDF-4, solo R; Python requiere el paquete netCDF4)"
if Rscript R/run_all.R "$TMP/nc4.yml" "$TMP/out_nc4_R" > "$TMP/nc4.R.log" 2>&1; then echo "  OK   R lee NetCDF-4"; else echo "  FALLA R con NetCDF-4"; status=1; fi
"$PY" python/exceedance.py --config "$TMP/nc4.yml" --output-dir "$TMP/out_nc4_Py" > "$TMP/nc4.py.log" 2>&1
grep -q "NetCDF-4/HDF5" "$TMP/nc4.py.log" && echo "  OK   Python sin netCDF4 da mensaje claro" || \
  { grep -q "Listo" "$TMP/nc4.py.log" && echo "  OK   Python lee NetCDF-4 (netCDF4 instalado)" || { echo "  FALLA Python con NetCDF-4"; status=1; }; }
echo "=== Resultado global: $([ $status -eq 0 ] && echo 'todo correcto' || echo 'HAY FALLAS')"
exit $status
