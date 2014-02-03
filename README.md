emulator
========

Netezza access:
In order to access Netezza on Rush for Phase 1, you need to load the custom built netezza module.  To do this, first log onto Rush and clone this repo into your home directory. Then copy the emulator/privatemodules directory to your home directoy.  Edit the netezza module file inside that directoy and add your database, username, and password.  Change the permissions on this file so that only you can read it and load the module.

$ cd ~
$ git clone https://github.com/kylemarcus/emulator.git
$ cp -r emulator/privatemodules .
* edit the netezza module file *
$ module load use.own
$ module load netezza

Phase 1:
The first phase of the emulator is to load the files into netezza and extract some information from the database. This is done by running the src/scripts/emulator_netezza.sh file.  Please look at the comments in this file to see exactly how to run it.  There will be some initial files that will have to be present first before running this script.

$ emulator/src/scripts/emulator_netezza.sh

Phase 2:
The second phase is running the MPI emulator code. You may need to change some of the variables defined in the makefile and slurm script file. The makefile is set up to use Intel compilers and the Intel MKL so these modules also need to be loaded.  Please note that you might also have to make some path changes in the DEFINE statements in the main cpp code before you compile.

$ module load intel
$ module load intel-mpi
$ module load mkl
$ cd emulator/src
$ make
$ sbatch emulator.slurm

Phase 3:
For phase 3 you must run src/scripts/merge.sh

$ emulator/src/scripts/merge.sh
