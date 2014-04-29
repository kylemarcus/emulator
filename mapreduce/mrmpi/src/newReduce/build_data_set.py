#!/usr/bin/env python

import sys
sys.path.append('.')
from neighbor import *
from parser import *
import re

# May 30 10:13 am - Message
# Need to rename this file as mapper1.py
# the purpose of this file is to following:
# 1.Read uncertain _input_list and store the contents as 'StructuredDataSet' class instance called samples. This can be encclosed in a function and
# there is currently no need to find sampl-sample neighbours.
# 2.Read uniq_coords as 'DataRead' and store it as 'StructuredDataSet' instance just like and samples and enclose in a function.
# 3.This is the real mapper part. If hadoop works right every node would have received a slice of the resamples file. Read the contents of 
# resamples file slice using 'DataRead' and then find its neighbors with samples. The samples function will return the samples instance 
# so find neighbors of resamples form among resamples and throw out the key value pairs. The output would have a column od resample keys against
# resample keys. That is the samples are keys and resamples are values.
# 4. repeat step 3 for phm. So call function 'uniq_coords_fun' to get the StructuredDataSet object and find neighbors of phm points
# from uniq_coords. Output would be uniq_coords keys against phm keys.
# Oh by the way the output should be stdout

# hadoop will disseminate data based on samples (keys). Reducer1.py will just write down all the resample keys for a sample to a unique file.



def Build_Sample_Object():
# Read uncertain input list as samples
    fp = open('/panasas/scratch/shivaswa/emulator/my_hadoop/uncertain_input_list.txt','r')
    skip_lines = int(re.search(r'\d+',fp.readline()).group())

    junk = [fp.readline() for i in range(skip_lines+2)]

    unc_inp_file = DataRead()
    sample_data = unc_inp_file.load_from_file(fp)
    samples = StructuredDataSet(transpose(sample_data))
    samples.scale()                             #Scaling the samples. Note that the original list is modified in the samples instance
    samples.hash_map()                          #To convert the samples list to key value form
    fp.close()
    return samples
#    sample_sample_neigh = samples.find_neighbors(samples.get_hashed_data(),0.2)       # sample-sample neighbor
#    disp_dictionary(sample_sample_neigh)


def Build_Uniq_Coords_Object():
    fp = open('/panasas/scratch/shivaswa/emulator/my_hadoop/Nx_Ny_group.txt','r')
    X_Y = re.findall('([\d.]+)',fp.readline())
    Nx_Ny_file = DataRead()
    grid = StructuredDataSet( transpose(Nx_Ny_file.load_from_file(fp)) )
    grid.hash_map(0)
    grid_dict = grid.get_hashed_data()
    fp.close()

    fp = open('/panasas/scratch/shivaswa/emulator/my_hadoop/uniq_coords.txt','r')
    uniq_coords_keys = [coords for coords in gen_X_Y_coord(fp,X_Y,grid_dict)]
    uniq_coords = StructuredDataSet(transpose(uniq_coords_keys))
#    uniq_coords.scale(0)
    uniq_coords.hash_map(0)
    fp.close()
    return uniq_coords


def gen_X_Y_coord(fp,X_Y,grid_dict):
    for key in fp:
        sample_group = 0
        match = re.search(r'(\d+)_(\d+)_(\d+)',key)
        my_group = int(match.group(1))
        if my_group != sample_group:   
            Nx = grid_dict[my_group][0]
            Ny = grid_dict[my_group][1]
            X = map( ( lambda num : (2*num + 0.5)*( float(X_Y[1]) - float(X_Y[0]) )/(2*Nx) + float(X_Y[0]) ), range(Nx) )
            Y = map( ( lambda num : (2*num + 0.5)*( float(X_Y[3]) - float(X_Y[2]) )/(2*Ny) + float(X_Y[2]) ), range(Ny) )
        row = int(match.group(2))-1
        col = int(match.group(3))-1
#       yield [match.group(0),X[col],Y[row]]
        yield [key.rstrip('\n'),X[col],Y[row]]
