#!/usr/bin/env bash
# Creates file for the Adaptive Rejection sampling function.

cat ./functions/adaptive_rejection_sampling_main.R > ars_all.R
cat ./functions/filter.R >> ars_all.R
cat ./functions/make_lowerupper_bound.R >> ars_all.R
cat ./functions/make_z.R >> ars_all.R
cat ./functions/sample_upper_bound.R >> ars_all.R
cat ./functions/update_x.R >> ars_all.R


