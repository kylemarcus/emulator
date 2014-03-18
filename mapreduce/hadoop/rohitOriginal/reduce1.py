#!/usr/bin/env python

import sys
sys.path.append('.')
import re
import operator as op
import itertools as it
from neighbor import *
from parser import gen_reverse

def gen_key_value(fp):
    for line in fp:
        match = re.findall('[\d._]+',line)
        yield [match[0], int(match[1])]

if __name__ == '__main__':
    fp = sys.stdin
    uniq_coords_phm_key = [key_value for key_value in gen_key_value(fp)]
    uniq_phm_dict = {}
    for uniq,full in it.groupby(sorted(uniq_coords_phm_key,key=op.itemgetter(0)),key=op.itemgetter(0)):
        phm_list = [res for res in gen_reverse(list(full))]
        uniq_phm_dict[uniq] = phm_list
#    disp_dictionary(uniq_phm_dict)
    
#    fp = open('/scratch/reduce1_output','w')
    for key in uniq_phm_dict:
	print key + '\t' + ' '.join( str(uniq_phm_dict[key]).rstrip('\]').lstrip('\[').split(',') )
#	fp.write(key + ' : ' + ' '.join( str(uniq_phm_dict[key]).rstrip('\]').lstrip('\[').split(',') ) + '\n' )

#    fp.close()
