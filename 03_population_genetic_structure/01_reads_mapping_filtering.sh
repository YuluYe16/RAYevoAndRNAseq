#!/usr/bin/env bash
set -euo pipefail

# Trim adapters and low-quality bases from paired-end reads.
trimmomatic PE -threads 32 -phred33 ray_1.fq.gz ray_2.fq.gz ray_1.paired.clean.fq.gz ray_1.unpaired.clean.fq.gz ray_2.paired.clean.fq.gz ray_2.unpaired.clean.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Align reads and retain properly paired primary alignments.
bwa mem -t 32 -M -v 1 -R '@RG\tID:ray\tSM:ray\tLB:ray' reference.genome.fa ray_1.paired.clean.fq.gz ray_2.paired.clean.fq.gz | samtools view -@ 4 -h -f 2 -F 2816 -o ray.reference.sam

# Identify and remove reads with at least 40 bp of cumulative soft clipping.
perl sam_get_soft_clip_readname.v2.pl ray.reference.sam 40 > ray.softclip.readnames.txt
perl sam_mask_with_read_name.pl ray.reference.sam ray.softclip.readnames.txt > ray.reference.clean.sam

# Sort alignments, remove PCR duplicates, and index the final BAM file.
samtools sort -@ 32 -O BAM ray.reference.clean.sam | samtools rmdup - ray.reference.clean.sorted.rmdup.bam
samtools index ray.reference.clean.sorted.rmdup.bam
