#!/bin/bash
#SBATCH --job-name=basCing
#SBATCH -o run_%j.out
#SBATCH --mail-type=ALL
#!SBATCH --mail-user=krishnap561997@gmail.com

#SBATCH --nodes=1
#SBATCH --ntasks=1                  # ONE task only
#SBATCH --cpus-per-task=32           # OpenMP threads (adjust!)
#!SBATCH --mem-per-cpu=4gb
#SBATCH -t 4-02:10:00

#SBATCH --account=bala1s
#SBATCH --qos=bala1s
#!SBATCH --partition=hpg-dev

# -------------------------
# Move to submit directory
# -------------------------
cd $SLURM_SUBMIT_DIR

date
hostname

# -------------------------
# Modules 
# -------------------------
module purge
module load gcc/14.2.0
module load imagemagick/7.0.8-20
module load hdf5

# -------------------------
# OpenMP environment
# -------------------------
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export KMP_STACKSIZE=1g

# Optional: debugging
# export OMP_DISPLAY_ENV=TRUE

# -------------------------
# Run
# -------------------------
echo "HPC_HDF5_LIB=$HPC_HDF5_LIB"
export LD_LIBRARY_PATH=$HPC_HDF5_LIB:$LD_LIBRARY_PATH

./$1
