#!/bin/bash
# Submit the no-fg-marg control on dr6plus_fiducial_baseline.
sbatch /home/jiaqu/dr6plus_lenslike/runs/fg_marg_test/cobaya.sh \
       /home/jiaqu/dr6plus_lenslike/runs/fg_marg_test/mcmc_control_baseline.yaml
