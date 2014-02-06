import os
import re
import glob
import argparse

"""
Creates output files needed for phase 1 of the emulator based on pileheight 
records and the uncertain input list file.

Outputs:
phData.txt      - this was rohit.xt
phMetaData.txt  - this was meta_data.txt

Inputs:
pileheights: /projects/gmfg/montserrat-2012/d_20cells/output_files_20cells
unvertain input list: /projects/gmfg/montserrat-2012/d_20cells/uncertain_input_list.txt
"""

parser = argparse.ArgumentParser()
parser.add_argument("uncertainInputList", help="Uncertain Input List file")
parser.add_argument("inputDir", help="pileheight input directoy")
parser.add_argument("outputDir", help="output directory")
args = parser.parse_args()

# uncertain input list matrix
u = []

# read in uncertain_input_list.txt
with open(args.uncertainInputList, 'r') as uilFp:
    
    # read header
    for i in range(6):
        samples = uilFp.readline()
    
    # read in initial condition matrix
    for i in range(int(samples)):
        u += uilFp.readline().split()

# write output files based on pileheight records
with open(os.path.join(args.outputDir, 'phData.txt'), 'w') as dataFp, \
     open(os.path.join(args.outputDir, 'phMetaData.txt'), 'w') as metaDataFp:
    
    # go through each pileheight record
    for phRecord in glob.iglob(os.path.join(args.inputDir, 'pileheight*')):
        
        filenumber = int(os.path.splitext(phRecord)[1][1:])
        print filenumber, phRecord
        
        with open(phRecord, 'r') as phRecordFp:

            # read header
            Nx, xstart, xend = re.findall(r'[\d.]+', phRecordFp.readline())
            Ny, ystart, yend = re.findall(r'[\d.]+', phRecordFp.readline())
            phRecordFp.readline()
            
            # read pileheight matrix
            pileHeight = []
            for i in range(int(Ny)):
                pileHeight += phRecordFp.readline().split()
            
            # write to data file
            for i in range(int(Ny)):
                for j in range(int(Nx)):
                    dataFp.write(str(filenumber) + ',' + \
                                 str(i+1) + ',' + \
                                 str(j+1) + ',' + \
                                 pileHeight[i*int(Nx)+j] + '\n')
            
            # write to meta data file
            i = (filenumber-1)*4
            metaDataFp.write(str(filenumber) + ',' + \
                             str(Nx) + ',' + \
                             str(Ny) + ',' + \
                             str(xstart) + ',' + \
                             str(xend) + ',' + \
                             str(ystart) + ',' + \
                             str(yend) + ',' + \
                             u[i] + ',' + \
                             u[i+1] + ',' + \
                             u[i+2] + ',' + \
                             u[i+3] + '\n')