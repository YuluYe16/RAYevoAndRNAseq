#!/usr/bin/env bash
set -euo pipefail
[[ $# -eq 3 ]] || { echo "usage: $0 BAM_LIST maf005|nomaf NEW_OUTPUT_DIR" >&2; exit 2; }
bam_list=$1
case "$2" in
  maf005) min_maf=0.05 ;;
  nomaf)  min_maf=-1 ;;
  *) echo "Choose maf005 or nomaf" >&2; exit 2 ;;
esac
[[ -r "$bam_list" ]]
awk 'NF {n++; if (++seen[$0]>1) dup++} END {exit(n!=724 || dup!=0)}' "$bam_list"
out=$3
mkdir -- "$out"
for chr in {1..12}; do
    angsd -uniqueOnly 1 -remove_bads 1 -minMapQ 30 -only_proper_pairs 1 \
        -nThreads 8 -b "$bam_list" -gl 2 -doMajorMinor 1 -doGlf 2 \
        -SNP_pval 1e-6 -doMaf 1 -minMaf "$min_maf" -minInd 580 \
        -doHaploCall 1 -doCounts 1 -minQ 13 -baq 0 -C 0 -trim 0 \
        -r "Chr${chr}:" -out "$out/asn724.Chr${chr}" \
        > "$out/asn724.Chr${chr}.log" 2>&1
done
