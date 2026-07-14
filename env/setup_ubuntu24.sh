#!/usr/bin/env bash
set -euo pipefail

# Bootstrap a reproducible PANTERA build environment on Ubuntu 24.04 / WSL2.
# The core executable needs MPI + PETSc. Python packages are for plotting/GUI
# post-processing scripts used by the case folders.

if ! command -v sudo >/dev/null 2>&1; then
   echo "sudo is required for apt package installation." >&2
   exit 1
fi

BASE_PACKAGES=(
   build-essential
   gfortran
   make
   pkg-config
   openmpi-bin
   libopenmpi-dev
   python3
   python3-tk
   python3-numpy
   python3-matplotlib
   python3-pandas
   python3-scipy
)

echo "==> Updating apt package index"
sudo apt-get update

echo "==> Installing compilers, MPI, PETSc, and Python post-processing packages"
if ! sudo apt-get install -y "${BASE_PACKAGES[@]}" libpetsc-real3.19-dev; then
   echo "==> libpetsc-real3.19-dev was unavailable; trying generic petsc-dev"
   sudo apt-get install -y "${BASE_PACKAGES[@]}" petsc-dev
fi

echo "==> Checking environment"
"$(dirname "$0")/check_env.sh"

echo "==> Building pantera.exe"
make -j"${JOBS:-$(nproc)}"

echo
echo "Environment is ready."
echo "Example:"
echo "  cp pantera.exe 2026/H2_low_recycling/"
echo "  cd 2026/H2_low_recycling"
echo "  mpirun -np 10 ./pantera.exe"
