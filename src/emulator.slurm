#!/bin/sh
#SBATCH --partition=general-compute
#SBATCH --time=15:00:00
#SBATCH --nodes=43
#SBATCH --ntasks-per-node=12
#SBATCH --mem=48000
#SBATCH --job-name="emulator"
#SBATCH --output=emulator.out
#SBATCH --mail-user=kmarcus2@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --exclusive

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR

cd $SLURM_SUBMIT_DIR
echo "working directory = "$SLURM_SUBMIT_DIR

module load intel
module load intel-mpi
module load mkl
module load use.own
module load netezza

ulimit -s unlimited

echo "Launch emulator with srun"

export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
srun -n 512 ./final_emulator_fork > /panasas/scratch/kmarcus2/montesrrat.out

echo done
