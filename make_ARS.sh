#!/usr/bin/env bash
# Creates file for the Adaptive Rejection sampling function.


cat adaptive_rejection_sampling_main.R >> ARS.R
cat filter.R >> ARS.R
cat make_lowerupper_bound.R >> ARS.R
cat make_z.R >> ARS.R
cat sample_upper_bound.R >> ARS.R
cat update_x.R >> ARS.R

