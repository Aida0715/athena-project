#!/usr/bin/env bash
#SBATCH --job-name=athena-smoke81
#SBATCH --partition=M-test-md
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=16G
#SBATCH --time=00:10:00
#SBATCH --output=cfca-smoke81-%j.out
#SBATCH --error=cfca-smoke81-%j.err

set -Eeuo pipefail

# ~/athena-project から、oneAPI環境を有効化した状態でsbatchする。
ROOT="${SLURM_SUBMIT_DIR}"
INPUT="${ROOT}/inputs/hydro/Toyouchi_CFCA/athinput.cfca-smoke-test81"
OUTDIR="/work/${USER}/results/cfca-smoke81-${SLURM_JOB_ID}"

test -x "${ROOT}/bin/athena"
test -r "${INPUT}"
mkdir -p "${OUTDIR}"

# 実際に使う実行ファイルとinputを保存する。
cp "${ROOT}/bin/athena" "${OUTDIR}/athena"
cp "${INPUT}" "${OUTDIR}/athinput"

{
  date -Is
  hostname
  git -C "${ROOT}" rev-parse HEAD
  git -C "${ROOT}" status --short
  module list
  sha256sum "${OUTDIR}/athena" "${OUTDIR}/athinput"
} > "${OUTDIR}/run-info.txt" 2>&1

export OMP_NUM_THREADS=1
cd "${OUTDIR}"

srun ./athena -i athinput > athena.log 2>&1
