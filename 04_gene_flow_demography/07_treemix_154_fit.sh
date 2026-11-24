#!/usr/bin/env bash
set -euo pipefail
[[ $# -eq 2 && -r "$1" ]] || { echo "usage: $0 COUNTS_GZ NEW_OUTPUT_DIR" >&2; exit 2; }
IFS= read -r header < <(gzip -dc -- "$1") || [[ -n "$header" ]]
[[ " $header " =~ [[:space:]]Outgroup[[:space:]] ]]
mkdir -- "$2"
for m in {0..10}; do
    for r in {1..10}; do
        prefix="$2/treemix.m${m}.r${r}"
        treemix -global -noss -seed "$((261113000+100*m+r))" -bootstrap 1000 \
            -k 2000 -i "$1" -m "$m" -root Outgroup -o "$prefix" \
            > "$prefix.log" 2> "$prefix.err"
    done
done
