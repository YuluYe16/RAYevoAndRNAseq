
# Build a de novo transposable-element library and annotate repeats.
EDTA.pl --genome ray.genome.fa --species Rice --overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --threads 32

# Soft-mask repeats with the EDTA-derived library.
RepeatMasker -e rmblast -s -gff -nolow -no_is -norna -xsmall -lib ray.genome.fa.mod.EDTA.TElib.fa -pa 20 ray.genome.fa

# Align paired-end RNA-seq reads to the genome.
hisat2-build -p 12 ray.genome.fa ray
hisat2 --dta -p 16 -x ray -1 ray_RNA_1.fq.gz -2 ray_RNA_2.fq.gz -S ray.rnaseq.sam
samtools sort -@ 16 -o ray.rnaseq.sorted.bam ray.rnaseq.sam
samtools index ray.rnaseq.sorted.bam

# Assemble genome-guided transcripts from RNA-seq alignments.
Trinity --max_memory 100G --genome_guided_bam ray.rnaseq.sorted.bam --genome_guided_max_intron 20000 --CPU 32 --output ray_trinity

# Predict protein-coding genes using RNA-seq and homologous-protein evidence.
braker.pl --genome=ray.genome.fa.masked --prot_seq=homologous_proteins.fa --bam=ray.rnaseq.sorted.bam --species=ray --threads=8 --workingdir=ray_braker --gff3 --nocleanup

# Align assembled transcripts to the masked genome with PASA.
Launch_PASA_pipeline.pl -c ray.alignAssembly.config --ALIGNERS gmap --MAX_INTRON_LENGTH 20000 -C -R -g ray.genome.fa.masked -t ray_trinity/Trinity-GG.fasta --CPU 32

# Load BRAKER gene models into the PASA database.
Load_Current_Gene_Annotations.dbi -c ray.alignAssembly.update.config -g ray.genome.fa.masked -P ray_braker/braker.gff3

# Update gene structures and UTRs using transcript evidence.
Launch_PASA_pipeline.pl -c ray.annotCompare.config -A -g ray.genome.fa.masked -t ray_trinity/Trinity-GG.fasta --CPU 32
