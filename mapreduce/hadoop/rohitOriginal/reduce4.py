#!/usr/bin/env python
import sys
sys.path.append('.')
import re



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
