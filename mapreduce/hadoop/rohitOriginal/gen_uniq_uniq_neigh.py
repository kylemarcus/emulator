import sys
from parser import *
from build_data_set import *

# This mapper first builds uniq_key object from uniq_coords.txt and then reads the uniq_coord_keys again from stdin to find 
# neighbors amog the uniq_coords i.e. this mapper finds neighbors of uniq_coords among uniq_coords.
# There is no reducer required for this. All generated files simply need to be concatenated.

if __name__ == '__main__':

    uniq_coords_obj = Build_Uniq_Coords_Object()
    uniq_coords_key_dict = uniq_coords_obj.get_hashed_data()

    fp = sys.stdin
#    count = 0;
    for key in fp:
#	count += 1
	my_key = key.rstrip('\n')
	neigh_list =  uniq_coords_obj.get_neighbors( uniq_coords_obj.get_parameter(my_key), 200 )
	my_str =  str(neigh_list).rstrip('\]').lstrip('\[')
#	print 'count = ' , count
	print my_key + ' : ' + ' '.join((re.split("\, |\'",my_str)))
	

