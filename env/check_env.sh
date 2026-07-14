#!/usr/bin/env bash
set -euo pipefail

missing=0

check_cmd() {
   local cmd="$1"
   if command -v "$cmd" >/dev/null 2>&1; then
      printf "OK   %-12s %s\n" "$cmd" "$(command -v "$cmd")"
   else
      printf "MISS %-12s\n" "$cmd"
      missing=1
   fi
}

check_cmd mpifort
check_cmd mpirun
check_cmd make
check_cmd pkg-config
check_cmd python3

if pkg-config --exists PETSc; then
   echo "OK   PETSc        $(pkg-config --modversion PETSc)"
   echo "     cflags:      $(pkg-config --cflags PETSc)"
   echo "     libs:        $(pkg-config --libs PETSc)"
else
   echo "MISS PETSc pkg-config entry"
   echo "     Install libpetsc-real3.19-dev or petsc-dev."
   missing=1
fi

python3 - <<'PY'
import importlib.util

mods = ["numpy", "matplotlib", "pandas", "tkinter"]
missing = [m for m in mods if importlib.util.find_spec(m) is None]
if missing:
    print("WARN optional Python post-processing modules missing:", ", ".join(missing))
else:
    print("OK   Python post-processing modules")
PY

if [ "$missing" -ne 0 ]; then
   echo
   echo "Environment check failed."
   exit 1
fi

echo
echo "Environment check passed."
