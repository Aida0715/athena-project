#!/usr/bin/env bash
#SBATCH --job-name=cfca-model71
#SBATCH --partition=M-bulk-md
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=cfca-model71-%j.out
#SBATCH --error=cfca-model71-%j.err

set -Eeuo pipefail

# oneAPI環境を有効化し、リポジトリ直下から投入する。
# 初回：sbatch run/CFCA/model71.sh
# 再開：sbatch run/CFCA/model71.sh /work/aidarj/results/前回の出力ディレクトリ
ROOT="${SLURM_SUBMIT_DIR}"
PREVIOUS="${1:-}"
WALL_TIME="${ATHENA_WALL_TIME:-11:45:00}"
OUTDIR="/work/${USER}/results/cfca-model71-${SLURM_JOB_ID}"

if [[ -z "${PREVIOUS}" ]]; then
  SOURCE="/work/${USER}/results/cfca-smoke81-335545"
  INPUT="${ROOT}/inputs/hydro/Toyouchi_CFCA/athinput.cfca-model71"
  test -x "${SOURCE}/athena"
  test -r "${INPUT}"
else
  SOURCE="$(cd -- "${PREVIOUS}" && pwd)"
  test -x "${SOURCE}/athena"
  test -s "${SOURCE}/Toyouchi.final.rst"
  test -r "${SOURCE}/athinput"
fi

mkdir "${OUTDIR}"
cp "${SOURCE}/athena" "${OUTDIR}/athena"
cp "${SOURCE}/run-info.txt" "${OUTDIR}/source-run-info.txt"
cp "${ROOT}/run/CFCA/model71.sh" "${OUTDIR}/job-script.sh"

if [[ -z "${PREVIOUS}" ]]; then
  cp "${INPUT}" "${OUTDIR}/athinput"
  args=(-i athinput -t "${WALL_TIME}")
else
  cp "${SOURCE}/athinput" "${OUTDIR}/athinput"
  cp "${SOURCE}/Toyouchi.final.rst" "${OUTDIR}/restart-source.rst"
  args=(-r restart-source.rst -t "${WALL_TIME}")
fi

export OMP_NUM_THREADS=1
unset DISPLAY
cd "${OUTDIR}"

{
  date -Is
  hostname
  echo "source_directory = ${SOURCE}"
  echo "previous_run = ${PREVIOUS}"
  echo "mpi_ranks = ${SLURM_NTASKS}"
  echo "repository_commit = $(git -C "${ROOT}" rev-parse HEAD)"
  module list
  sha256sum athena athinput
  if [[ -n "${PREVIOUS}" ]]; then
    sha256sum restart-source.rst
  fi
  printf 'command = srun --cpu-bind=cores ./athena'
  printf ' %q' "${args[@]}"
  printf '\n'
} > run-info.txt 2>&1

set +e
srun --cpu-bind=cores ./athena "${args[@]}" > athena.log 2>&1
status=$?
set -e

{
  echo "end_time = $(date -Is)"
  echo "exit_code = ${status}"
} >> run-info.txt

if (( status != 0 )); then
  exit "${status}"
fi

test -s Toyouchi.final.rst
sha256sum Toyouchi.final.rst >> run-info.txt
