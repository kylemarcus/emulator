#!/usr/bin/env python
import sys
sys.path.append('.')
from neighbor import *
from parser import print_key_value
from build_data_set import *


if __name__ == '__main__':
    filename = sys.argv[1]

    fp = open('/user/shivaswa/my_hadoop/montserrat_take2_vol_dir_bed_int.phm','r')
    max_min_dict = get_max_min(fp)
    fp.close()

# Build phm data 
    fp = sys.stdin
    phm_read_file = DataRead()
    phm_data = phm_read_file.load_from_file(fp)
    tran_phm_data = transpose(phm_data)
#    scaled_phm_data = [sub_list for sub_list in gen_scaled_data(tran_phm_data,max_min_dict)]
    phm = StructuredDataSet(tran_phm_data)
    phm.hash_map(0)
    phm_dict = phm.get_hashed_data()
    disp_dictionary(phm_dict)

# Build Uniq coords data
    Uniq_Coords = Build_Uniq_Coords_Object()
    uniq_coords_dict = Uniq_Coords.get_hashed_data()
#    disp_dictionary(uniq_coords_dict)
    phm_uniq_coords_neigh_dict = phm.find_neighbors(uniq_coords_dict,200)

# Write the count of uniq_id neighbors for each phm_id to a file whose name is the hostname of the cpu it is running on
    filename = filename+"_"+str(phm_dict.keys()[0])
    fp = open(filename,"w")
    fp.write(str(phm_dict.keys()[0]))
#    for key in phm_uniq_coords_neigh_dict:
#	my_str = str(key) + " " + str(len(phm_uniq_coords_neigh_dict[key])) + "\n"
#	fp.write(my_str)
    fp.close()

    for key in phm_uniq_coords_neigh_dict:
        print_key_value(key,phm_uniq_coords_neigh_dict[key])
