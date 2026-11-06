#!/usr/bin/env bash
#SBATCH --job-name=athena-restart81
#SBATCH --partition=M-test-md
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=16G
#SBATCH --time=00:10:00
#SBATCH --output=cfca-restart81-%j.out
#SBATCH --error=cfca-restart81-%j.err

set -Eeuo pipefail

# oneAPI環境を有効化してからsbatchする。
SOURCE="/work/${USER}/results/cfca-smoke81-335545"
OUTDIR="/work/${USER}/results/cfca-restart81-${SLURM_JOB_ID}"

test -x "${SOURCE}/athena"
test -r "${SOURCE}/Toyouchi.final.rst"
mkdir -p "${OUTDIR}"

cp "${SOURCE}/athena" "${OUTDIR}/athena"
cp "${SOURCE}/Toyouchi.final.rst" "${OUTDIR}/restart-source.rst"
cp "${SOURCE}/run-info.txt" "${OUTDIR}/source-run-info.txt"

{
  date -Is
  hostname
  echo "source_run = ${SOURCE}"
  echo "restart_override = time/nlim=20"
  module list
  sha256sum "${OUTDIR}/athena" "${OUTDIR}/restart-source.rst"
} > "${OUTDIR}/run-info.txt" 2>&1

export OMP_NUM_THREADS=1
cd "${OUTDIR}"
srun ./athena -r restart-source.rst time/nlim=20 > athena.log 2>&1
