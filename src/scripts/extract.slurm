#!/bin/sh
#SBATCH --partition=general-compute
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name="extract"
#SBATCH --mail-user=kmarcus2@buffalo.edu
#SBATCH --mail-type=ALL

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR

cd $SLURM_SUBMIT_DIR
echo "working directory = "$SLURM_SUBMIT_DIR

module load python-epd/7.1.2

ulimit -s unlimited

echo "Running python extract.py"

python extract.py /projects/gmfg/montserrat-2012/d_20cells/output_files_20cells /panasas/scratch/kmarcus2 &>/panasas/scratch/kmarcus2/extract.eo

echo done
