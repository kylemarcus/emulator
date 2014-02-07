CC=mpiicpc
FLAGS=-I$(MKLROOT)/include
LIBS=-L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm

installDir = /panasas/scratch/kmarcus2/emulator/run
slurmScript = emulator.slurm

objects = final_emulator_fork

all: $(objects)

$(objects): *.cpp
	$(CC) $(FLAGS) $@.cpp -o $@ $(LIBS)

install:
	mkdir -p $(installDir)
	mv $(objects) $(installDir)
	cp $(slurmScript) $(installDir)

clean:
	rm $(objects)
