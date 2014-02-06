import os
import re
import glob
import argparse

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

with open(os.path.join(args.inputDir, 'rohit.txt'), 'w') as outputFp:
    
    for phRecord in glob.iglob(os.path.join(args.inputDir, 'pileheight*')):
        
        pileHeight = []

        filenumber = int(os.path.splitext(phRecord)[1][1:])
        print phRecord, filenumber
        
        with open(phRecord,'r') as phRecordFp:

            # read header
            Nx, xstart, xend = re.findall(r'[\d.]+', phRecordFp.readline())
            Ny, ystart, yend = re.findall(r'[\d.]+', phRecordFp.readline())
            fp.readline()
            
            # read pileheight matrix
            for i in range(int(Ny)):
                pileHeight += phRecordFp.readline().split()
            
            # write to output file
            for i in range(int(Ny)):
                for j in range(int(Nx)):
                    outputFp.write( str(filenumber) + ',' + \
                                    str(i+1) + ',' + \
                                    str(j+1) + ',' + \
                                    str(pileHeight[i*Nx+j]) + '\n' )
                    