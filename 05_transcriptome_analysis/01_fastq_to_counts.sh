
# Build the HISAT2 index for the reference genome.
hisat2-build -p 8 ray.genome.fa ray_hisat2_index

# Align paired-end RNA-seq reads and sort the alignments.
hisat2 -p 8 --dta --new-summary -x ray_hisat2_index -1 ray_R1.fq.gz -2 ray_R2.fq.gz 2> ray.hisat2.log | samtools sort -@ 4 -o ray.bam -

# Index the BAM file and summarize mapping statistics.
samtools index ray.bam
samtools flagstat ray.bam > ray.flagstat.txt

# Count unstranded paired-end fragments over exon features grouped by Parent.
featureCounts -T 16 -p --countReadPairs -s 0 -F GFF -t exon -g Parent -a ray.gene.gff3 -o ray.gene_counts.txt ray.bam
