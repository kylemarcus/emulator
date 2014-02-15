CC=mpiicpc
FLAGS=-g -I$(BOOST_INC_DIR) -I$(MKLROOT)/include
LIBS=-L$(BOOST_LIB_DIR) -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -lboost_program_options

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
