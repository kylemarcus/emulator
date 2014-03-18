#!/usr/bin/env python
import sys
sys.path.append('.')
import re
import operator as op
import itertools as it
from neighbor import *
from parser import gen_key_value
from parser import gen_reverse


# In this file sam stands for sample and res stands for resample

if __name__ == '__main__':
    fp = sys.stdin
    sam_res = [key_value for key_value in gen_key_value(fp)]
    sam_res_dict = {}
    for uniq,full in it.groupby(sorted(sam_res,key=op.itemgetter(0)),key=op.itemgetter(0)):
	res_list = [res for res in gen_reverse(list(full))]
	sam_res_dict[uniq] = res_list

    for key in sam_res_dict:
        print key , ' '.join( str(sam_res_dict[key]).rstrip('\]').lstrip('\[').split(',') )

