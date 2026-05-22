#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=20
#SBATCH --output=/scratch/jiaqu/mpi_output_%j.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jq247@cam.ac.uk

cd $SLURM_SUBMIT_DIR
export DISABLE_MPI=false
export PYTHONPATH="/home/jiaqu/dr6plus_lenslike:$PYTHONPATH"
export PYTHONPATH=/home/jiaqu/cobaya:$PYTHONPATH

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
export COBAYA_USE_FILE_LOCKING=False

YAML=${1:?usage: sbatch cobaya.sh <yaml_path>}

# Resolve absolute output dir from the YAML's `output:` field and clean stale locks
OUT_DIR=$(python3 -c "
import sys, yaml
with open('$YAML') as f:
    cfg = yaml.safe_load(f)
out = cfg.get('output', '')
import os
print(os.path.dirname(out) if out else '')
")
if [[ -n "$OUT_DIR" && -d "$OUT_DIR" ]]; then
    echo "Cleaning lock files in $OUT_DIR"
    find -L "$OUT_DIR" -name "*lock*" -delete 2>/dev/null
fi

srun -n 4 cobaya-run --force "$YAML"
