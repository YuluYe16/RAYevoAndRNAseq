#!/usr/bin/env bash
set -euo pipefail
[[ $# -eq 4 && -d "$1" && -r "$2" && -r "$3" ]] || {
    echo "usage: $0 CHR_VCF_DIR SETS_TSV TREE_NEWICK NEW_OUTPUT_DIR" >&2; exit 2;
}
for c in {1..12}; do [[ -r "$1/Chr${c}.vcf.gz" ]]; done
mkdir -- "$4"
prefixes=()
for c in {1..12}; do
    prefix="$4/Chr${c}"
    Dsuite Dtrios -k 5 -t "$3" -o "$prefix" "$1/Chr${c}.vcf.gz" "$2" \
        > "$prefix.log" 2>&1
    prefixes+=("$prefix")
done
Dsuite DtriosCombine -t "$3" -o "$4/fb" "${prefixes[@]}" > "$4/combine.log" 2>&1
Dsuite Fbranch "$3" "$4/fb_combined_tree.txt" > "$4/fbranch.tsv"
