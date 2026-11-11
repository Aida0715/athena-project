#!/usr/bin/env bash
set -Eeuo pipefail
# Delete ONLY the abandoned model71 long-run experiment, keeping smoke/benchmark data.
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
command -v squeue >/dev/null
command -v scancel >/dev/null
queue="$(squeue --noheader --user="${USER}" --format='%i %j')"
jobs=()
while read -r id name; do
  if [[ "${name}" == cfca-model71 && "${id}" =~ ^[0-9]+$ ]]; then
    jobs+=("${id}")
  fi
done <<< "${queue}"
if ((${#jobs[@]})); then
  scancel "${jobs[@]}"
  # Wait for Slurm to stop writers before removing their directories.
  for attempt in {1..30}; do
    remaining="$(squeue --noheader --jobs="$(IFS=,; echo "${jobs[*]}")" --format='%i')"
    [[ -z "${remaining}" ]] && break
    sleep 2
  done
  if [[ -n "${remaining}" ]]; then
    echo 'ERROR: model71 jobs still stopping; rerun cleanup later.' >&2
    exit 1
  fi
fi
python3 - "${ROOT}" "/work/${USER}/results" <<'PY'
from pathlib import Path
import re
import shutil
import sys
root, results = map(Path, sys.argv[1:])
for folder, pattern in ((results, r'cfca-model71-\d+'),
                        (root, r'cfca-model71-\d+\.(out|err)')):
    if not folder.is_dir():
        continue
    for p in sorted(folder.iterdir()):
        if not re.fullmatch(pattern, p.name):
            continue
        if p.is_dir() and not p.is_symlink():
            shutil.rmtree(p)
        else:
            p.unlink()
        print('Deleted:', p)
PY
