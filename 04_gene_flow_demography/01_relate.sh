
# Prepare phased haplotypes, samples, ancestral states, and population labels for Relate.
for chromosome in {1..12}; do
    PrepareInputFiles.sh --haps "Chr${chromosome}.relate.haps" --sample "Chr${chromosome}.relate.sample" --ancestor ancestral.fa -o "Chr${chromosome}_input" --poplabels ray.poplabels
done

# Infer chromosome-specific genealogies using the recombination maps.
for chromosome in {1..12}; do
    Relate --mode All -m 6.5e-09 -N 30000 --haps "Chr${chromosome}_input.haps.gz" --sample "Chr${chromosome}_input.sample.gz" --map "Chr${chromosome}.map" --annot "Chr${chromosome}_input.annot" --seed 1 -o "ray_chr${chromosome}"
done

# Re-estimate branch lengths and population sizes across all chromosomes.
EstimatePopulationSize.sh -i ray -m 6.5e-09 --poplabels ray.poplabels --years_per_gen 1 --first_chr 1 --last_chr 12 --num_iter 5 --bins 2,7,0.1 --threads 32 --seed 1 -o ray_reest

# Plot population-size and relative cross-coalescence-rate trajectories.
Rscript plotPopSize.R -i ray_reest.pairwise.coal -o ray --min 100 --max 1000000 --years_per_gen 1 --focal_pop RAY

# Convert re-estimated Relate outputs to tskit tree sequences.
for chromosome in {1..12}; do
    RelateFileFormats --mode ConvertToTreeSequence -i "ray_reest_chr${chromosome}" -o "Chr${chromosome}"
done
