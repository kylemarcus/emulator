#!/bin/bash

# Workflow script for mr-mpi
# Kyle Marcus
# Summer 14

case "$1" in

    install)

        echo installing MR-MPI

        module load intel
        module load intel-mpi
        
        cd ~
        wget http://www.cs.sandia.gov/~sjplimp/tars/mrmpi.tar.gz
        tar xf mrmpi.tar.gz
        cd mrmpi-*/src
        make -f Makefile.shlib mpicc

        echo installing PyPar

        module load python-epd/7.1.2

        cd ~
        git clone https://github.com/daleroberts/pypar.git
        cd pypar
        python setup.py bulid
        python setup.py install --user

        echo installing mr-mpi custom modules
        
        cd ~
        git clone https://github.com/kylemarcus/emulator.git 
        mkdir -p privatemodules/mrmpi-python
        cp emulator/privatemodule/mrmpi-python/22Nov13 privatemodules/mrmpi-python

        ;;

    run)

        echo run ...
        ;;

    *)

        echo $"Usage: $0 {install|run}"
        exit 1

esac
