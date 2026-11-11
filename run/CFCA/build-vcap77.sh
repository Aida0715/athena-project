#!/usr/bin/env bash
#SBATCH --job-name=cfca-vcap77-build
#SBATCH --partition=M-test-md
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --mem=16G
#SBATCH --time=00:30:00
set -Eeuo pipefail
BUNDLE="${1:?Usage: build-vcap77.sh BUNDLE}"
cd "${BUNDLE}/source"
{
  date -Is
  hostname
  module list
  cat "${BUNDLE}/source-commit.txt"
  icpx --version
  mpiicpx --version
} > "${BUNDLE}/build-info.txt" 2>&1
JOBS="${SLURM_CPUS_PER_TASK}" bash scripts/build/cfca.sh
cp bin/athena "${BUNDLE}/athena"
srun --ntasks=1 "${BUNDLE}/athena" -c >> "${BUNDLE}/build-info.txt" 2>&1
(cd "${BUNDLE}" && sha256sum athena > executable.sha256)
