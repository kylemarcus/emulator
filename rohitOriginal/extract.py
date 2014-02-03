import numpy as np
import os
import fnmatch as fn
import re

# creates the file rohit.txt, this file contains a list of the header information of
# all the pileheight records

# pileheights: /projects/gmfg/montserrat-2012/d_20cells/output_files_20cells

# rohit.txt contains:
# number Nx X[0] X[1] Ny Y[0] Y[1]

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

#	for i in range(Ny):
#		line = fp.readline()
#		h = line.split()
#		H = H + h
	
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

#	print len(H)
#	del(H)
#	fp.close()


#	if index == 1:
#		fp = open('Check_data.txt','w')
#	else:
#		fp = open('Check_data.txt','a')

	
#	for i in range(6):
#		if i==0:
#			p = str(filenumber) + ',' + 'Nx' + ',' + str(data[i])
#		elif i==1:
#			p = str(filenumber) + ',' + 'Xstart' + ',' + str(data[i])
#		elif i==2:
#			p = str(filenumber) + ',' + 'Xend' + ',' + str(data[i])
#		elif i==3:
#			p = str(filenumber) + ',' + 'Ny' + ',' + str(data[i])
#		elif i==4:
#			p = str(filenumber) + ',' + 'Ystart' + ',' + str(data[i])
#		elif i==5:
#			p = str(filenumber) + ',' + 'Yend' + ',' + str(data[i])

				
#		p = str(filenumber) + ' ' + str(data[0]) + ' ' + str(data[1])+ ' ' + str(data[2])+ ' ' + str(data[3])+ ' ' + str(data[4])+ ' ' + str#(data[5])

#		fp.write(p)
		fp.write('\n')
	
	fp.close()
