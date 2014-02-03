import numpy as np
import os
import fnmatch as fn
import re

# this file creates the meta_data.txt
Number = []
data = []
index = 0

# uncertain_input_list is the initial conditions that are being looked at,
# each col in the file represents a different attribute
fp = open('uncertain_input_list.txt','r');

# read header, skips to the 6th line
for i in range(6):
	samples = fp.readline()

# read the number of samples specified on the 6th line
Nsamples = int(samples)

U = []

# read the data matrix into U line by line
for i in range(Nsamples):
	samples = fp.readline()
	samples = samples.split()
	U = U + samples

i = 0
files = fn.filter(os.listdir(os.getcwd()),'pileheight*')

# go through each pileheight record
for F in files:
        index = index + 1
        H = []
        data = []

        file = str(F)
        print file
        filenumber = int(re.search(r'[\d]+',file).group())
        print filenumber
        fp = open(file,'r')
        line = fp.readline()

        xline = re.findall(r'[\d.]+',line)
        line = fp.readline()
        yline = re.findall(r'[\d.]+',line)
        data = data + xline + yline
        Nx = int(xline[0])
        Ny = int(yline[0])
        xstart = float(xline[1])
        ystart = float(yline[1])
        xend = float(xline[2])
        yend = float(yline[2])
	fp.close()

	if index == 1:
                fp = open('meta_data.txt','w')
        else:
                fp = open('meta_data.txt','a')

	i = (filenumber-1)*4;
	fp.write(str(filenumber)+','+str(Nx)+','+str(Ny)+','+str(xstart)+','+str(xend)+','+str(ystart)+','+str(yend))
	fp.write(','+str(U[i])+','+str(U[i+1])+','+str(U[i+2])+','+str(U[i+3]))
	fp.write('\n')
	fp.close()
