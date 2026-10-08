#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT}"

JOBS="${JOBS:-4}"
if [[ ! "${JOBS}" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: JOBS must be a positive integer" >&2
  exit 2
fi

python3 configure.py \
  --prob=Toyouchi \
  --coord=cartesian \
  --eos=isothermal \
  --flux=hlld \
  --grav=mg \
  --nghost=4 \
  -b \
  -mpi \
  --cxx=g++ \
  --mpiccmd=mpicxx

make clean
make -j "${JOBS}"
