#!/usr/bin/env python
import sys
sys.path.append('.')
import re


def gen_output(LIMIT,my_list):
  fp = sys.stdin
  for line in fp:
    match = re.findall( '([\d_]+)', line )
    key = match[0]
#    iter = 0
    try:
	id = my_list.index(key)
    except:
	continue

    for iter,phm in enumerate(match):
       if iter == 0:
	continue
       if iter <= LIMIT:
	 print phm," ",key
	 iter += 1
       else:
	 break


if __name__ == '__main__':
# 2 arguments are required for running this script : 1.Sample number   2. Limit on the neighbors
# This file just prints out the phm_id and uniq_id against it. THe sample number is required so because each downsampled data
# has different uniq_ids. THis reads the respective downsampled files and generates output for only those uniq_ids that exist in the downsampled
# file. The LIMIT parameter is to keep a limit on the number of neighbors considered, just for computational time consideration. 
# Right now I use a limit of 35 (took cue from Netezza results.)

  downsample_number = int(sys.argv[1])
  filename = 'down_sample_files/downsample%06d.txt'%downsample_number
  LIMIT = int( sys.argv[2])
  uniq_id_list = []
  fp = open(filename,"r")
  for line in fp:
     match = re.findall( '([\d_.]+)', line )
     uniq_id_list = uniq_id_list + [match[1]]

  gen_output(LIMIT,uniq_id_list)
