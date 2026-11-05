#!/usr/bin/env bash
set -Eeuo pipefail

# 実行前に同じターミナルで開発環境を有効化する：
# source /work/opt/local/bin/enable-oneapi.sh

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT}"

JOBS="${JOBS:-4}"
if [[ ! "${JOBS}" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: JOBS must be a positive integer" >&2
  exit 2
fi

if ! command -v icpx >/dev/null ||
   ! command -v mpiicpx >/dev/null; then
  echo "ERROR: Activate the CFCA Intel oneAPI environment first." >&2
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
  --cxx=icpx \
  --mpiccmd=mpiicpx

make clean
make -j "${JOBS}"
