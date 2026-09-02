
# Estimate alleles, genotype likelihoods, haploid genotypes, and allele frequencies at filtered SNP sites.
for chromosome in {1..12}; do
    angsd -uniqueOnly 1 -remove_bads 1 -minMapQ 30 -only_proper_pairs 1 -doHaploCall 1 -nThreads 32 -doCounts 1 -bam ray.bamlist -ref reference.genome.fa -out "ray.Chr${chromosome}" -r "Chr${chromosome}" -sites "Chr${chromosome}.clean.site" -doMajorMinor 1 -GL 1 -doMaf 1
done
