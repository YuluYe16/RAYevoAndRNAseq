
# Detect terminal telomeric-repeat arrays.
seqtk telo -m CCCTAAA ray.genome.fa > ray.telomere.txt

# Assess conserved single-copy ortholog completeness.
busco -i ray.genome.fa -l embryophyta_odb10 -o ray_busco -m genome --cpu 16 --offline

# Map Illumina reads and assess short-read coverage.
bwa index ray.genome.fa
bwa mem -t 32 ray.genome.fa ray_clean_1.fq.gz ray_clean_2.fq.gz | samtools sort -@ 8 -o ray.illumina.sorted.bam
samtools index ray.illumina.sorted.bam
qualimap bamqc -bam ray.illumina.sorted.bam -nt 8 --java-mem-size=80G -c -outdir ray_illumina_qualimap

# Map HiFi reads and assess long-read coverage.
minimap2 -t 32 -ax map-hifi ray.genome.fa ray.hifi.fastq | samtools sort -@ 32 -o ray.hifi.sorted.bam
samtools index ray.hifi.sorted.bam
qualimap bamqc -bam ray.hifi.sorted.bam -nt 32 --java-mem-size=100G -c -outdir ray_hifi_qualimap

# Build an Illumina k-mer database and estimate consensus QV.
meryl k=19 count output ray_1.meryl ray_clean_1.fq.gz threads=12
meryl k=19 count output ray_2.meryl ray_clean_2.fq.gz threads=12
meryl union-sum output ray.meryl ray_1.meryl ray_2.meryl threads=32
merqury.sh ray.meryl ray.genome.fa ray_QV

# Identify candidate LTR retrotransposons for LAI estimation.
gt suffixerator -db ray.genome.fa -indexname ray -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index ray -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > ray.harvest.scn
LTR_FINDER_parallel -seq ray.genome.fa -threads 32 -harvest_out -size 1000000
cat ray.harvest.scn ray.genome.fa.finder.combine.scn > ray.rawLTR.scn
LTR_retriever -genome ray.genome.fa -inharvest ray.rawLTR.scn -threads 32

# Calculate the LTR Assembly Index.
LAI -genome ray.genome.fa -intact ray.genome.fa.pass.list -all ray.genome.fa.out -t 32
