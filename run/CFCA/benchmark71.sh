#!/usr/bin/env bash
#SBATCH --job-name=athena-bench71
#SBATCH --partition=M-test-md
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=32G
#SBATCH --time=00:15:00
#SBATCH --output=cfca-bench71-%j.out
#SBATCH --error=cfca-bench71-%j.err

set -Eeuo pipefail

# oneAPI環境を有効化し、リポジトリ直下から
# sbatch --ntasks=<並列数> run/CFCA/benchmark71.sh で投入する。
ROOT="${SLURM_SUBMIT_DIR}"
SOURCE="/work/${USER}/results/cfca-smoke81-335545"
INPUT="${ROOT}/inputs/hydro/Toyouchi_CFCA/athinput.cfca-benchmark71"
OUTDIR="/work/${USER}/results/cfca-bench71-${SLURM_NTASKS}r-${SLURM_JOB_ID}"

test -x "${SOURCE}/athena"
test -r "${INPUT}"
mkdir -p "${OUTDIR}"
cp "${SOURCE}/athena" "${OUTDIR}/athena"
cp "${INPUT}" "${OUTDIR}/athinput"
cp "${SOURCE}/run-info.txt" "${OUTDIR}/source-run-info.txt"
cp "${ROOT}/run/CFCA/benchmark71.sh" "${OUTDIR}/job-script.sh"

{
  date -Is
  hostname
  echo "mpi_ranks = ${SLURM_NTASKS}"
  echo "executable_source = ${SOURCE}/athena"
  git -C "${ROOT}" rev-parse HEAD
  module list
  sha256sum "${OUTDIR}/athena" "${OUTDIR}/athinput"
} > "${OUTDIR}/run-info.txt" 2>&1

export OMP_NUM_THREADS=1
unset DISPLAY
cd "${OUTDIR}"

start_seconds=$(date +%s)
set +e
srun --cpu-bind=cores ./athena -i athinput > athena.log 2>&1
status=$?
set -e

{
  echo "end_time = $(date -Is)"
  echo "elapsed_seconds = $(($(date +%s) - start_seconds))"
  echo "exit_code = ${status}"
} >> run-info.txt

exit "${status}"
