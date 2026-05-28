#!/bin/bash
# Submit the A_fg-marg chain (A_fg ~ U[-5, 5]) on dr6plus_fiducial_extended.
sbatch /home/jiaqu/dr6plus_lenslike/runs/fg_marg_test/cobaya.sh \
       /home/jiaqu/dr6plus_lenslike/runs/fg_marg_test/mcmc_afg_extended_u5.yaml
