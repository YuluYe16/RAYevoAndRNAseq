
# Calculate D statistics for all population trios in parallel.
DtriosParallel -k 200 --cores 64 -n K200 -t tree.nwk populations.txt ray.vcf.gz

# Calculate f-branch statistics on the population tree.
Dsuite Fbranch tree.nwk DTparallel_ray_K200_combined_tree.txt > ray_K200_tree.fbranch.txt

# Plot the f-branch matrix and population tree.
python3 dtools.py ray_K200_tree.fbranch.txt tree.nwk
