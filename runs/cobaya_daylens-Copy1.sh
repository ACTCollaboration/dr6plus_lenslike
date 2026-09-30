#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=20
#SBATCH --output=/scratch/iabril/chains/output/mpi_output_%j.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ia404@cam.ac.uk
#SBATCH --account=rrg-rbond-ac

cd $SLURM_SUBMIT_DIR
# export PYTHONPATH="/home/jiaqu/dr6plus_lenslike:$PYTHONPATH"

module load StdEnv/2023
module load aocl-lapack/5.1
module load openblas
module load gsl
module load openmpi
module load fftw
module load cfitsio
module load python
module load StdEnv/2023  intel/2023.2.1  openmpi/4.1.5
module load mpi4py/3.1.4

source /scratch/iabril/.virtualenv/likelihood/bin/activate

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export COBAYA_USE_FILE_LOCKING=False

export DISABLE_MPI=false

srun -n 1 cobaya-run daytime_lcdm_daylens_Llow.yaml

