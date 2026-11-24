#!/usr/bin/env bash
set -euo pipefail
[[ $# -eq 6 ]] || { echo "usage: $0 CHROM QUERY REF SAMPLE_MAP MAP NEW_OUTPUT_DIR" >&2; exit 2; }
[[ "$1" =~ ^Chr([1-9]|1[0-2])$ ]]
for f in "$2" "$3" "$4" "$5"; do [[ -r "$f" ]]; done
mkdir -- "$6"
rfmix -f "$2" -r "$3" -m "$4" -g "$5" --chromosome="$1" \
    -o "$6/$1" --n-threads=4 --random-seed=20261103 \
    -G 8 -s 100 -c 20 -t 100 -e 0 > "$6/$1.log" 2>&1
