#!/usr/bin/env bash
set -euo pipefail
[[ $# -eq 2 && -r "$1" ]] || { echo "usage: $0 SNP_FASTA NEW_OUTPUT_DIR" >&2; exit 2; }
awk '/^>/ {n++; if (++ids[$1]>1) dup++} END {exit(n!=552 || dup!=0)}' "$1"
mkdir -- "$2"
FastTree -nt -gtr "$1" > "$2/asn552.RAY14.K9Q07.gtr.fasttree.nwk" 2> "$2/FastTree.log"
