#!/usr/bin/env bash
set -Eeuo pipefail
# Run on CFCA after activating oneAPI. No build occurs on the login node.
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT}"
command -v sbatch >/dev/null
command -v icpx >/dev/null
command -v mpiicpx >/dev/null
# git archive includes committed sources only. Refuse uncommitted build changes.
if [[ -n "$(git status --porcelain -- src configure.py Makefile.in scripts/build/cfca.sh)" ]]; then
  echo 'ERROR: Commit build/source changes first.' >&2
  exit 2
fi
COMMIT="$(git rev-parse HEAD)"
BASE="/work/${USER}/results"
mkdir -p "${BASE}"
BUNDLE="$(mktemp -d "${BASE}/cfca-vcap77-$(date +%Y%m%d-%H%M%S)-XXXXXX")"
mkdir "${BUNDLE}/source" "${BUNDLE}/inputs" "${BUNDLE}/logs" "${BUNDLE}/results"
git archive "${COMMIT}" | tar -x -C "${BUNDLE}/source"
printf '%s\n' "${COMMIT}" > "${BUNDLE}/source-commit.txt"
for limit in 100 300 1000; do
  cp "inputs/hydro/Toyouchi_CFCA/athinput.cfca-vcap77-${limit}" "${BUNDLE}/inputs/"
done
cp run/CFCA/{build-vcap77,vcap77,submit-vcap77}.sh "${BUNDLE}/"
(cd "${BUNDLE}" && sha256sum inputs/* *.sh > setup.sha256)
echo "Bundle: ${BUNDLE}"
build_reply="$(sbatch --parsable --output="${BUNDLE}/logs/build-%j.out" \
  --error="${BUNDLE}/logs/build-%j.err" "${BUNDLE}/build-vcap77.sh" "${BUNDLE}")"
BUILD_ID="${build_reply%%;*}"
[[ "${BUILD_ID}" =~ ^[0-9]+$ ]]
printf 'build_job=%s\n' "${BUILD_ID}" | tee "${BUNDLE}/submission.txt"
array_reply="$(sbatch --parsable --dependency="afterok:${BUILD_ID}" \
  --kill-on-invalid-dep=yes --output="${BUNDLE}/logs/run-%A_%a.out" \
  --error="${BUNDLE}/logs/run-%A_%a.err" "${BUNDLE}/vcap77.sh" "${BUNDLE}")"
ARRAY_ID="${array_reply%%;*}"
[[ "${ARRAY_ID}" =~ ^[0-9]+$ ]]
printf 'array_job=%s\nindex_0=100 km/s\nindex_1=300 km/s\nindex_2=1000 km/s\n' \
  "${ARRAY_ID}" | tee -a "${BUNDLE}/submission.txt"
echo "Check: squeue -j ${BUILD_ID},${ARRAY_ID}"
