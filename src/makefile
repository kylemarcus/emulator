CC=mpiicpc
FLAGS=-I$(MKLROOT)/include
LIBS=-L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm

objects = final_emulator_fork

all: $(objects)

$(objects): *.cpp
	$(CC) $(FLAGS) $@.cpp -o $@ $(LIBS)

clean:
	rm $(objects)
