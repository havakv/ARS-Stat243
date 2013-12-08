#!/usr/bin/env bash
# Creates file for the Adaptive Rejection sampling function.



cat ./functions/adaptive_rejection_sampling_main.R > ARS.R
cat ./functions/filter.R >> ARS.R
cat ./functions/make_lowerupper_bound.R >> ARS.R
cat ./functions/make_z.R >> ARS.R
cat ./functions/sample_upper_bound.R >> ARS.R
cat ./functions/update_x.R >> ARS.R

