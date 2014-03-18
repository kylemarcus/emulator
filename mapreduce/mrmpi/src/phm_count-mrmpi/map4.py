#!/usr/bin/env python
import sys
sys.path.append('.')
import re

DOWN_SAMPLE_FILES='/panasas/scratch/kmarcus2/emulator/mrmpi/down_sample_files/downsample%06d.txt'
UNIQ_PHM='/user/kmarcus2/emulator/mapreduce/mrmpi/src/phm_count-mrmpi/uniq_phm.txt'

def gen_output(LIMIT,my_list):
  
  fp = open(UNIQ_PHM,"r")
  
  return_list = []
  
  for line in fp:
    
    match = re.findall( '([\d_]+)', line )
    key = match[0]
    
    try:
      id = my_list.index(key)
    except:
      continue

    for iter,phm in enumerate(match):
      if iter == 0:
        continue
      if iter <= LIMIT:
        return_list += [(phm,key)]
        iter += 1
      else:
        break
  
  return return_list

def main(filename, limit):
  
  fp = open(filename,"r")
  
  uniq_id_list = []
  
  for line in fp:
    match = re.findall( '([\d_.]+)', line )
    uniq_id_list += [match[1]]

  return gen_output(limit,uniq_id_list)


"""
Two arguments are required for running this script:
 1. Sample number
 2. Limit on the neighbors
 
This file just prints out the phm_id and uniq_id against it. THe sample number is required so because each downsampled data
has different uniq_ids. This reads the respective downsampled files and generates output for only those uniq_ids that exist in the downsampled
file. The LIMIT parameter is to keep a limit on the number of neighbors considered, just for computational time consideration. 
Right now I use a limit of 35 (took cue from Netezza results.)
"""
if __name__ == '__main__':
  
  downsample_number = int(sys.argv[1])
  filename = DOWN_SAMPLE_FILES%downsample_number
  limit = int(sys.argv[2])
  
  l = main(filename, limit)
  
  for t in l:
    print t[0], t[1]
  