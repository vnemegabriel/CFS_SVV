#!/bin/bash

# Slurm job options (job-name, compute nodes, job time)
#SBATCH --job-name=12-a12
#SBATCH --time=00:20:0
#SBATCH --nodes=1
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1

# Replace [budget code] below with your budget code (e.g. t01)
#SBATCH --account=e01-ICL-lah
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --distribution=block:block 
#SBATCH --hint=nomultithread 
#SBATCH --reservation=shortqos

# Setup the job environment (this module needs to be loaded before any other modules)
module load epcc-job-env

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically 
#   using threading.
export OMP_NUM_THREADS=1

export NEK_DIR=/work/e01/e01/mlahooti/nektarSpencer/build
export NEK_BUILD=$NEK_DIR/dist/bin
export LD_LIBRARY_PATH=/opt/gcc/9.3.0/snos/lib64:$NEK_DIR/ThirdParty/dist/lib:$NEK_DIR/dist/lib64:$LD_LIBRARY_PATH

# Launch the parallel job

srun $NEK_BUILD/IncNavierStokesSolver FP2D.xml session_ptub.xml -i Hdf5 &> runlog
