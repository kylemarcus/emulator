#!/usr/bin/env python
import sys
sys.path.append('.')
import re
import operator as op
import itertools as it
from neighbor import *
from parser import gen_key_value
from parser import gen_reverse

def gen_non_empty_list(my_dict):
    for key in my_dict:
	if len(my_dict[key]):
	    yield my_dict[key]
	    yield key

if __name__ == '__main__':

    sam_res_dist_dict = {}
    for i in range(2048):
	sam_res_dist_dict[i] = []


    fp = sys.stdin
    for line in fp:
	match = re.findall('[\d.]+',line)
	sample = int(match[0])
	resample = int(match[1])
	dist = float(match[2])
	sam_res_dist_dict[sample].append([resample,dist])

    my_iter = gen_non_empty_list(sam_res_dist_dict)
    for res_dist_list in my_iter:
	sample_num = my_iter.next()
	sorted_res_dist_list = sorted(res_dist_list, key = op.itemgetter(1))
	limit = 1200 if len(sorted_res_dist_list) > 1200 else len(sorted_res_dist_list)
	new_sam_res_dist = sorted_res_dist_list[0:limit]
	resample_tuple =  zip(*new_sam_res_dist)[0]
	print str(sample_num) + ' ' + ' '.join( re.split('\,',str(resample_tuple).lstrip('\(').rstrip('\)') ) )
	
#	    print str(sample_num) + ' ' + str(my_list[0]) + ' ' +  str(my_list[1])
