
# Run ADMIXTURE with cross-validation; repeat this command across K values and random seeds.
admixture --cv -j4 -s 12345 ray.bed 2 > ray.K2.log 2>&1
