#!/usr/bin/env bash
# Download the 53 participants (EEG, events) and their ICA weights from OpenNeuro
# (ds006979, version 1.1.5 or later) into data/ds006979: about 37 GB of EEG and
# 5.5 GB of ICA weights. Each file is compared with its size on OpenNeuro: complete
# files are skipped and incomplete ones (an interrupted download) downloaded again,
# so the script can simply be rerun after an interruption.
#
# If the vwm-erp-julia project already downloaded some participants, point both
# projects to the same folder instead of downloading again:
#   DS006979=/path/to/data/ds006979 bash scripts/download_data.sh
# If the ICA weights are already on disk in another folder (searched in all its
# subfolders), set ICA_WEIGHTS to that folder: they are then not downloaded, and the
# Julia scripts read them from there.
#   ICA_WEIGHTS=/path/to/ds006979_derivatives bash scripts/download_data.sh
#
# Usage (from the project folder):
#   bash scripts/download_data.sh

cd "$(dirname "$0")/.." || exit 1
ROOT=${DS006979:-data/ds006979}
BASE=https://s3.amazonaws.com/openneuro.org/ds006979
DER=derivatives/ReplicationPeterson2015
FAILED=0

# get URL OUT: download unless OUT already has the size of the file on OpenNeuro.
get() {
  local url=$1 out=$2 remote local_size
  remote=$(curl -sSIL --retry 3 "$url" | tr -d '\r' | awk 'tolower($1)=="content-length:"{n=$2} END{print n}')
  local_size=$( [ -f "$out" ] && wc -c < "$out" | tr -d ' ' || echo 0)
  [ -n "$remote" ] && [ "$local_size" = "$remote" ] && return 0
  if curl -sSfL --retry 3 -o "$out.part" "$url" && mv "$out.part" "$out"; then
    echo "  downloaded $(basename "$out")"
  else
    rm -f "$out.part"; echo "  failed: $out (rerun the script)"; FAILED=$((FAILED + 1))
  fi
}

for s in 01 02 03 04 07 08 09 10 11 12 13 14 15 16 17 18 19 21 22 23 24 25 28 30 \
         32 33 34 35 36 39 40 41 43 46 47 48 49 50 51 52 53 54 58 59 61 62 64 \
         69 70 72 73 74 76; do
  d=$ROOT/sub-$s/eeg; mkdir -p "$d"
  for f in eeg.vhdr eeg.vmrk eeg.eeg events.tsv; do
    get "$BASE/sub-$s/eeg/sub-${s}_task-visualworkingmemorytask_$f" \
        "$d/sub-${s}_task-visualworkingmemorytask_$f"
  done
  n=$((10#$s))
  # ICA weights already on disk elsewhere (ICA_WEIGHTS=/path/to/folder): not downloaded.
  if [ -n "${ICA_WEIGHTS:-}" ] && \
     [ -n "$(find "$ICA_WEIGHTS" -name "${n}_allERPs_optimized_weights.set" -size +1000k 2>/dev/null | head -1)" ]; then
    echo "sub-$s done (ICA weights found in ICA_WEIGHTS)"; continue
  fi
  w=$ROOT/$DER/sub-$s; mkdir -p "$w"
  get "$BASE/$DER/sub-$s/${n}_allERPs_optimized_weights.set" "$w/${n}_allERPs_optimized_weights.set"
  echo "sub-$s done"
done

echo "Participants: $(ls -d "$ROOT"/sub-* | wc -l | tr -d ' ') (expected 53)"
[ "$FAILED" -eq 0 ] && echo "All files complete." || echo "$FAILED file(s) failed: run the script again."
