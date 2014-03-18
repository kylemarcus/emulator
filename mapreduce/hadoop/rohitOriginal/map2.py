#!/usr/bin/env python
import sys
sys.path.append('.')
from build_data_set import *
from parser import *
import time

if __name__ == '__main__':

    Samples = Build_Sample_Object()
    samples_dict = Samples.get_hashed_data()
#    disp_dictionary(Samples.get_hashed_data())
#    Uniq_Coords = Build_Uniq_Coords_Object()
#    disp_dictionary(Uniq_Coords.get_hashed_data())


    fp = open('/user/shivaswa/my_hadoop/resamples.txt','r')
#    fp.readline()
    max_min_dict = get_max_min(fp)
    fp.close()

    fp = sys.stdin
    read_resamples_file = DataRead()
    resamples_data = read_resamples_file.load_from_file(fp)
    tran_resamples_data = transpose(resamples_data)
# Remember that resamples data has a column of weights too that needs to be read from the file

    prune_tran_resamples_data = tran_resamples_data[0:read_resamples_file.cols()-1]
#    scaled_resamples_data = [scale(sub_list, max_min_dict[key+1][0], max_min_dict[key+1][1]) for key,sub_list in enumerate(prune_tran_resamples_data)]

    scaled_resamples_data = [sub_list for sub_list in gen_scaled_data(prune_tran_resamples_data,max_min_dict)]

    Resamples = StructuredDataSet(scaled_resamples_data)
    Resamples.hash_map(0)
    resamples_dict = Resamples.get_hashed_data()
    disp_dictionary(resamples_dict)
#    resample_sample_neigh_dict = Resamples.find_neighbors(samples_dict,.12)

#    for key in resample_sample_neigh_dict:
#	print_key_value(key,resample_sample_neigh_dict[key])

# I wrote this neighbor generator function ( gen_find_neighbors_with_dist in neighbor.py ) specifically for this. It generates sample key followed 
# by resample
# key and the distance of the sample-resample set of parameters. I did this to extract only the top 1200 or so. because with a parametric scaled
# distance of 0.2 samples appear to have more than 14000 resamples as neighbors. Now unless someone wants to kill oneself with that sort
# of ridiculous number of calclations (and think of the cost of reduce operation at the end of emulator ) one would really want to use this.
# The idea is to filter and retain only top 1200 in the reduce2.py. I got the number 1200 by cheating....yes, snooping at the results obtained
# from netezza. The max number of resample neighbors for a sample there was ~1200. So I'm gong to fix 1200 as the limit. Oh boy 14000 - you must
# be kidding me. Once again the output is : Sample_key Resample_key Distance

    for result in Resamples.gen_find_neighbors_with_dist(samples_dict,.12):
	print ' '.join(re.split('\,',(str(result).lstrip('\[').rstrip('\]')) ) )
