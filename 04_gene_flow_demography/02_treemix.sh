
# Calculate population-stratified allele frequencies from the filtered VCF.
plink --vcf ray.vcf.gz --freq gz --within treemix.cluster --double-id --out ray --allow-extra-chr --set-missing-var-ids @:# --keep-allele-order

# Convert PLINK frequencies to TreeMix input format.
python2.7 plink2treemix.py ray.frq.strat.gz ray.treemix.frq.gz

# Repeat TreeMix analyses for zero to fifteen migration edges.
mkdir -p treemix_rep
for migration_edges in {0..15}; do
    for replicate in {1..10}; do
        treemix -global -noss -seed "$RANDOM" -bootstrap 1000 -k 2000 -i ray.treemix.frq.gz -m "$migration_edges" -root outgroup -o "treemix_rep/ray.m${migration_edges}.r${replicate}" > "treemix_rep/ray.m${migration_edges}.r${replicate}.log" 2> "treemix_rep/ray.m${migration_edges}.r${replicate}.err"
    done
done
