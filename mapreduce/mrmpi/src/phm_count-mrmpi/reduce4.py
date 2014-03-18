#!/usr/bin/env python
import sys
sys.path.append('.')
import re

def main(map_list):
  
  resample_count_dict = {}
  return_list = []
  
  for t in map_list:
    
    line = str(t[0]) + ' ' + str(t[1])
    
    match = re.findall( '([\d_]+)', line )
    key = match[0]
    
    if key in resample_count_dict:
      resample_count_dict[key] += 1
    else:
      resample_count_dict[key] = 1

  for key in resample_count_dict:
    return_list += [(key,resample_count_dict[key])]
  
  return return_list

if __name__ == '__main__':
  
  resample_count_dict = {}
  fp = sys.stdin
  for line in fp:
    match = re.findall( '([\d_]+)', line )
    key = match[0]
    if key in resample_count_dict:
      resample_count_dict[key] += 1
    else:
      resample_count_dict[key] = 1

  for key in resample_count_dict:
    print key," ",resample_count_dict[key]
