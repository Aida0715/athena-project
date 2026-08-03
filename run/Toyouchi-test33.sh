#!/usr/bin/env bash

set -Eeuo pipefail

# commit: cd2fe2c //本シミュレーションに対応するToyouchi.cppをgitの履歴から追跡可
# Toyouchi+23の設定に磁場Bz=1μGを入れてテスト

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ATHENA="${ATHENA_DIR:-$(cd -- "${SCRIPT_DIR}/.." && pwd)}"
WORK="${WORK_DIR:-/work/beta/aida}"
OUTDIR="${OUTPUT_DIR:-${WORK}/results/Toyouchi-test33}"  #毎回変更
INPUT="${INPUT_FILE:-${ATHENA}/inputs/hydro/Toyouchi_test/athinput.Toyouchi_25}"  #毎回変更
NPROC="${NPROC:-8}"  #MPI並列計算のコア数

if [[ ! "${NPROC}" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: NPROC must be a positive integer: ${NPROC}" >&2
  exit 2
fi

if [[ ! -x "${ATHENA}/bin/athena" ]]; then
  echo "ERROR: MPI executable not found or not executable: ${ATHENA}/bin/athena" >&2
  exit 2
fi

if [[ ! -r "${INPUT}" ]]; then
  echo "ERROR: input file is not readable: ${INPUT}" >&2
  exit 2
fi

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

RUN_NAME="$(basename "$0" .sh)"
RUN_LOG="${RUN_LOG:-${OUTDIR}/${RUN_NAME}.$(date +%Y%m%d-%H%M%S).log}"

# 実行（ファイル名はinputに任せる）
# 環境作成中のテスト計算用入力ファイルはToyouchi_testディレクトリに格納
# 並列計算のコア数は-np32のように指定する
unset DISPLAY  #GUIを使わない

{
  echo "start_time = $(date --iso-8601=seconds)"
  echo "executable = ${ATHENA}/bin/athena"
  echo "input      = ${INPUT}"
  echo "output_dir = ${OUTDIR}"
  echo "mpi_ranks  = ${NPROC}"
  echo "log_file   = ${RUN_LOG}"
} | tee -a "${RUN_LOG}"

# "$@" に Athena++ の上書き引数（例: time/nlim=10）を渡せる。
# PIPESTATUS[0] を使い、tee ではなく mpirun 本体の終了コードを返す。
set +e
mpirun --use-hwthread-cpus -np "${NPROC}" "${ATHENA}/bin/athena" \
  -i "${INPUT}" "$@" 2>&1 | tee -a "${RUN_LOG}"
status=${PIPESTATUS[0]}
set -e

echo "end_time   = $(date --iso-8601=seconds)" | tee -a "${RUN_LOG}"
echo "exit_code  = ${status}" | tee -a "${RUN_LOG}"

if (( status != 0 )); then
  echo "ERROR: Athena++/MPI terminated abnormally. See ${RUN_LOG}" >&2
fi

exit "${status}"

  #Problem generator:            Toyouchi
  #Coordinate system:            cartesian
  #Equation of state:            isothermal
  #Riemann solver:               hlle
  #Magnetic fields:              ON
  #Number of scalars:            0
  #Number of chemical species:   0
  #Special relativity:           OFF
  #General relativity:           OFF
  #Radiative Transfer:           OFF
  #Implicit Radiation:           OFF
  #Cosmic Ray Transport:         OFF
  #Cosmic Ray Diffusion:         OFF
  #Frame transformations:        OFF
  #Self-Gravity:                 multigrid
  #Super-Time-Stepping:          OFF
  #Chemistry:                    OFF
  #KIDA rates:                   OFF
  #ChemRadiation:                OFF
  #chem_ode_solver:              OFF
  #Debug flags:                  OFF
  #Code coverage flags:          OFF
  #Linker flags:                  
  #Floating-point precision:     double
  #Number of ghost cells:        4
  #MPI parallelism:              ON
  #OpenMP parallelism:           OFF
  #FFT:                          OFF
  #HDF5 output:                  OFF
  #Compiler:                     g++
  #Compilation command:          g++  -O3 -std=c++11

