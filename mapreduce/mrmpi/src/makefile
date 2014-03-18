all: newEmulator

newEmulator: matrixop.h new_build_emulator.h newEmulator.cpp
	icpc -Wl,-rpath=$(MKL)/lib/intel64,--start-group -I$(MKL)/include -L$(MKL)/lib/intel64 -lmkl_sequential -lmkl_intel_lp64 -lmkl_core -Wl,--end-group -lpthread newEmulator.cpp -o newEmulator
