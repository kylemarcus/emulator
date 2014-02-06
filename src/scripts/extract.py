import os
import re
import argparse
import fnmatch as fn
import numpy as np

"""
Creates the file rohit.txt, which has information on pileheight records

pileheights: /projects/gmfg/montserrat-2012/d_20cells/output_files_20cells

rohit.txt contains:
 sample col row pileheight
"""

parser = argparse.ArgumentParser()
parser.add_argument("inputDir", help="pileheight input directoy")
parser.add_argument("outputDir", help="rohit.txt output directory")
args = parser.parse_args()

Number = []
data = []
index = 0
files = fn.filter(os.listdir(os.getcwd()),'pileheight*')

for F in files:
	index = index + 1
	H = []
	data = []

	file = str(F)
	print file
	filenumber = int(re.search(r'[\d]+',file).group()) #gets the ending timestamp number of filename
	print filenumber
	fp = open(file,'r')
	line = fp.readline()

	xline = re.findall(r'[\d.]+',line) #gets first 3 numbers into array
	line = fp.readline()
	yline = re.findall(r'[\d.]+',line) #gets the next 3 numbers into array
	data = data + xline + yline
	Nx = int(xline[0])
	Ny = int(yline[0])
	xstart = float(xline[1])
	ystart = float(yline[1])
	xend = float(xline[2])
	yend = float(yline[2])
	fp.readline()
	
	fp.close()

	Number = Number + [filenumber]

	if index == 1:
		fp = open('rohit.txt','w')
	else:
		fp = open('rohit.txt','a')


	for i in range(Ny):
		for j in range(Nx):
			p = str(filenumber) + ',' + str(i+1) + ',' + str(j+1) + ',' + str(H[i*Nx + j])
			fp.write(p)
			fp.write('\n')
		fp.write('\n')
	
	fp.close()
