# Align the reference and ray genomes for structural-variation analysis.
minimap2 -ax asm5 -t 10 --eqx -o reference_vs_ray.sam reference.genome.fa ray.genome.fa

# Identify syntenic regions and structural rearrangements.
syri -c reference_vs_ray.sam -r reference.genome.fa -q ray.genome.fa --no-chrmatch -F S --prefix reference_vs_ray_

# Plot synteny and rearrangements between the two genomes.
plotsr --sr reference_vs_ray_syri.out --genomes genomes.txt --chrord chrord.txt --itx --nodup --notr -H 8.5 -W 10 -f 8 -o reference_vs_ray.synteny.pdf
