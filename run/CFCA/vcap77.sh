#!/usr/bin/env bash
#SBATCH --job-name=cfca-vcap77
#SBATCH --partition=M-bulk-md
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-2
set -Eeuo pipefail
# Submit using submit-vcap77.sh, which freezes inputs and builds one shared executable.
BUNDLE="${1:?Usage: vcap77.sh BUNDLE [PREVIOUS_RESULT]}"
PREVIOUS="${2:-}"
limits=(100 300 1000)
LIMIT="${limits[${SLURM_ARRAY_TASK_ID:?}]}"
OUTDIR="${BUNDLE}/results/v${LIMIT}-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir "${OUTDIR}"
(cd "${BUNDLE}" && sha256sum -c executable.sha256)
cp "${BUNDLE}/athena" "${OUTDIR}/athena"
cp "${BUNDLE}/inputs/athinput.cfca-vcap77-${LIMIT}" "${OUTDIR}/athinput"
cp "${BUNDLE}/vcap77.sh" "${OUTDIR}/job-script.sh"
cp "${BUNDLE}/build-info.txt" "${OUTDIR}/build-info.txt"
args=(-i athinput -t 11:45:00)
if [[ -n "${PREVIOUS}" ]]; then
  # A single array index must be selected for restart, e.g. --array=1 for 300 km/s.
  cmp "${OUTDIR}/athinput" "${PREVIOUS}/athinput"
  cmp "${OUTDIR}/athena" "${PREVIOUS}/athena"
  test -s "${PREVIOUS}/Toyouchi.final.rst"
  cp "${PREVIOUS}/Toyouchi.final.rst" "${OUTDIR}/restart-source.rst"
  args=(-r restart-source.rst -t 11:45:00)
fi
export OMP_NUM_THREADS=1
unset DISPLAY
cd "${OUTDIR}"
{
  date -Is
  hostname
  cat "${BUNDLE}/source-commit.txt"
  echo "velocity_ceiling_kms = ${LIMIT}"
  echo "previous_result = ${PREVIOUS}"
  echo "mpi_ranks = ${SLURM_NTASKS}"
  module list
  sha256sum athena athinput
  if [[ -n "${PREVIOUS}" ]]; then sha256sum restart-source.rst; fi
  printf 'command = srun --cpu-bind=cores ./athena'
  printf ' %q' "${args[@]}"
  printf '\n'
} > run-info.txt 2>&1
set +e
srun --cpu-bind=cores ./athena "${args[@]}" > athena.log 2>&1
status=$?
set -e
printf 'end_time = %s\nexit_code = %s\n' "$(date -Is)" "${status}" >> run-info.txt
if (( status != 0 )); then exit "${status}"; fi
test -s Toyouchi.final.rst
sha256sum Toyouchi.final.rst >> run-info.txt
# COMPLETED also includes a clean wall-time stop; record whether 20 kyr was reached.
if grep -q 'Terminating on time limit' athena.log; then
  echo 'target_status = physical time limit reached' >> run-info.txt
else
  echo 'target_status = inspect log; restart may be required' >> run-info.txt
fi
