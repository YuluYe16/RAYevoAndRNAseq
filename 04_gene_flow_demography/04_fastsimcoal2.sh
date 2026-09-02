
# Retain biallelic SNPs from callable neutral regions on the twelve chromosomes.
bcftools view -r Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8,Chr9,Chr10,Chr11,Chr12 -S samples.keep.txt -R neutral_callable.bed -m2 -M2 -v snps -Oz -o ray.neutral.vcf.gz ray.vcf.gz
tabix -p vcf ray.neutral.vcf.gz

# Preview feasible multidimensional SFS projection sizes.
python3 easySFS.py -i ray.neutral.vcf.gz -p ray.sample_pop.tsv -a -f --preview

# Generate a folded multidimensional SFS using the selected projections.
python3 easySFS.py -i ray.neutral.vcf.gz -p ray.sample_pop.tsv -o ray_easySFS --prefix ray_proj8 -a -f --proj 8,8,8,8,8,8

# Run one fastsimcoal2 optimization; repeat this command in separate replicate directories.
fsc28 --cores 8 --numBatches 12 -t ray.tpl -e ray.est --msfs --multiSFS -M -q -n 100000 -L 50 -s 0

# Select the replicate with the highest estimated likelihood.
python3 fsc_bestrun.py -i fsc_results/ray -o best_runs/ray

# Calculate AIC from the best-likelihood output and matching model definition.
Rscript calculateAIC.R ray
