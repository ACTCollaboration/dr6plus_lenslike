#!/bin/bash
# Submit the afg_extended_bao chain (lensing + DESI DR2 BAO) on dr6plus_fiducial_extended.
sbatch /home/jiaqu/dr6plus_lenslike/runs/fg_marg_test/cobaya.sh \
       /home/jiaqu/dr6plus_lenslike/runs/fg_marg_test/mcmc_afg_extended_bao.yaml
