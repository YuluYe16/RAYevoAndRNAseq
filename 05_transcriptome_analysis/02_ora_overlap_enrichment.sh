
# Identify genes whose bodies overlap OrA-derived genomic intervals.
bedtools intersect -u -a ray.genes.bed -b ray.OrA_tracts.bed > ray.genes_OrA_overlap.bed

# Identify genes whose bodies plus 2-kb flanks overlap OrA-derived intervals.
bedtools intersect -u -a ray.genes_flank2k.bed -b ray.OrA_tracts.bed > ray.genes_flank2k_OrA_overlap.bed

# Calculate overlap fold enrichment, one-sided hypergeometric P values, and BH-adjusted FDR.
Rscript ora_enrichment_statistics.R ray.OrA_overlap_input.tsv ray.OrA_enrichment.tsv
