import os
import itertools

def grouper(n, iterable):
    args = [iter(iterable)] * n
    return ([e for e in t if e != None] for t in itertools.izip_longest(*args))

inDir = "/panasas/scratch/kmarcus2/emulator/my_hadoop/down_sample_files/"
#inDir = "/panasas/scratch/kmarcus2/emulator/my_hadoop/down_sample_files_2000/"
groups = 6

files = os.listdir(inDir)
sorted_files = sorted(files, key=lambda x: int(x[10:16]))

grouped_files = list(grouper(groups, sorted_files))

for group in grouped_files:
    
    start = int(group[0][10:16])
    end = int(group[len(group)-1][10:16])

    out =  "#!/bin/sh\n"
    out += "#SBATCH --partition=general-compute\n"
    out += "#SBATCH --time=15:00:00\n"
    out += "#SBATCH --nodes=1\n"
    out += "#SBATCH --ntasks-per-node=12\n"
    out += "#SBATCH --mem=48000\n"
    out += "#SBATCH --job-name=\"emulator-"+str(start)+".."+str(end)+"\"\n"
    out += "#SBATCH --output=\"emulator-"+str(start)+".."+str(end)+".out\"\n"
    out += "#SBATCH --mail-user=kmarcus2@buffalo.edu\n"
    out += "#SBATCH --mail-type=ALL\n"
    out += "#SBATCH --exclusive\n"
    out += "\n"
    out += "echo \"SLURM_JOBID=\"$SLURM_JOBID\n"
    out += "echo \"SLURM_JOB_NODELIST\"=$SLURM_JOB_NODELIST\n"
    out += "echo \"SLURM_NNODES\"=$SLURM_NNODES\n"
    out += "echo \"SLURMTMPDIR=\"$SLURMTMPDIR\n"
    out += "\n"
    out += "cd /scratch\n"
    out += "rm file*\n\n"
    #out += "cd $SLURM_SUBMIT_DIR\n"
    #out += "echo \"working directory = \"$SLURM_SUBMIT_DIR\n"
    out += "cd /user/kmarcus2/emulator/mapreduce/mrmpi/src\n"
    out += "\n"
    out += "module load intel\n"
    out += "module load python-epd/7.1.2\n"
    out += "\n"
    out += "ulimit -s unlimited\n"
    out += "\n"
    out += "echo \"Launch emulator "
    for f in group:
        out += str(int(f[10:16])) + " "
    out += "\"\n\n"
    for f in group:
        out += "python emulator.py " + inDir + f + " | ./newEmulator &\n"
    out += "\n"
    out += "wait\n"
    out += "echo done with emulator\n"
    out += "echo transfering output files\n"
    out += "cd /scratch\n"
    #out += "cp /scratch/file* /panasas/scratch/kmarcus2/emulator/my_hadoop/emulator_output_using_2000\n"
    out += "tar cfz - file* | tar xfz - -C /panasas/scratch/kmarcus2/emulator/my_hadoop/emulator_output\n"
    out += "rm file*\n"
    out += "\n"
    out += "echo done"

    fp = open("emulator-"+str(start)+".."+str(end)+"-slurm","w")
    fp.write(out)
    fp.close()
    
