#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=20
#SBATCH --output=/scratch/r/rbond/jiaqu/mpi_output_%j.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jq247@cam.ac.uk
#SBATCH --account=rrg-rbond-ac

cd $SLURM_SUBMIT_DIR
export DISABLE_MPI=false

module load StdEnv/2023
module load aocl-lapack/5.1
module load openblas
module load gsl
module load openmpi
module load fftw
module load cfitsio
module load python
source /home/jiaqu/.bashrc
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export COBAYA_USE_FILE_LOCKING=false
export PYTHONPATH="/home/jiaqu/dr6plus_lenslike:$PYTHONPATH"


srun -n 4 cobaya-run act_lcdm.yaml

