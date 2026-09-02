# Convert PacBio HiFi BAM reads to FASTQ.
bam2fastq -u -o ray.hifi ray.bam

# Assemble HiFi reads and identify plant telomeric repeats.
hifiasm -t 32 --telo-m CCCTAAA -o ray ray.hifi.fastq > ray.hifiasm.log 2>&1

# Convert the primary-contig GFA file to FASTA.
awk '/^S/{print ">"$2; print $3}' ray.bp.p_ctg.gfa > ray.bp.p_ctg.fa

# Summarize contig-level assembly statistics.
assembly-stats ray.bp.p_ctg.fa > ray.bp.p_ctg.stat

# Scaffold primary contigs against the IR64-T2T reference genome.
ragtag.py scaffold IR64-T2T.fasta ray.bp.p_ctg.fa -u -t 32 --aligner minimap2 -f 20000 --remove-small -o ray_primary
