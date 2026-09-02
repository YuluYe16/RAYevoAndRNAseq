# RAYevoAndRNAseq

Concise command-line workflows used for genome evolution and RNA-seq analyses in RAY rice.

## Contents

- `01_genome_assembly_annotation`: HiFi genome assembly, quality assessment, repeat annotation, and gene annotation.
- `02_genome_alignment_phylogeny`: whole-genome alignment, synteny visualization, orthology inference, and species-tree reconstruction.
- `03_population_genetic_structure`: read mapping and filtering, ANGSD SNP analysis, PCA, and ADMIXTURE.
- `04_gene_flow_demography`: Relate, TreeMix, Dsuite F-branch, and fastsimcoal2 workflows.
- `05_transcriptome_analysis`: RNA-seq alignment, read counting, and OrA-tract overlap enrichment.

## Usage

The scripts provide the principal commands and parameters used in the analyses. The example prefix `ray` represents a generic sample or dataset and should be replaced as needed.

Input data, reference genomes, population labels, demographic model files, software-specific configuration files, and custom helper scripts are not included. Users should prepare these files according to the corresponding software documentation and their own data structure.

The scripts are intended as concise methodological records rather than a single end-to-end executable pipeline.
