# Set MIN_IND to half of the population sample count; set depth limits to 0.5 and 2.5 times the expected total depth.

# Estimate genotype likelihoods and the folded site-frequency spectrum for ray.
angsd -bam ray.bamlist -ref reference.genome.fa -anc reference.genome.fa -out ray -nThreads 10 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -C 50 -baq 1 -minMapQ 30 -minQ 30 -minInd "${MIN_IND}" -setMinDepth "${MIN_DEPTH}" -setMaxDepth "${MAX_DEPTH}" -doCounts 1 -GL 1 -doSaf 1

# Generate the corresponding SAF file for the comparison population.
angsd -bam comparison.bamlist -ref reference.genome.fa -anc reference.genome.fa -out comparison -nThreads 10 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -C 50 -baq 1 -minMapQ 30 -minQ 30 -minInd "${COMPARISON_MIN_IND}" -setMinDepth "${COMPARISON_MIN_DEPTH}" -setMaxDepth "${COMPARISON_MAX_DEPTH}" -doCounts 1 -GL 1 -doSaf 1

# Estimate the folded one-dimensional SFS.
realSFS ray.saf.idx -P 10 -fold 1 > ray.sfs

# Calculate nucleotide-diversity statistics from the folded SFS.
realSFS saf2theta ray.saf.idx -sfs ray.sfs -fold 1 -outname ray
thetaStat do_stat ray.thetas.idx -win 100000 -step 100000 -outnames ray.thetastat.win100k
thetaStat do_stat ray.thetas.idx -win 10000 -step 10000 -outnames ray.thetastat.win10k
thetaStat do_stat ray.thetas.idx -win 2000 -step 2000 -outnames ray.thetastat.win2k

# Estimate the folded two-dimensional SFS for ray and a comparison population.
realSFS ray.saf.idx comparison.saf.idx -P 10 -fold 1 > ray_comparison.ml

# Calculate global and 100-kb windowed pairwise FST.
realSFS fst index ray.saf.idx comparison.saf.idx -sfs ray_comparison.ml -fold 1 -whichFst 1 -fstout ray_comparison
realSFS fst stats ray_comparison.fst.idx > ray_comparison.global_fst.txt
realSFS fst stats ray_comparison.fst.idx -win 100000 -step 100000 > ray_comparison.win100k_fst.txt
