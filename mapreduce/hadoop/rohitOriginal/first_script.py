from parser import *

# This script should be the first one to run (hence the name) and can be run on the front end. Make sure to run this on a single node.
# This script performs two operations. First is to read the meta data from all pileheight records and group the samples based on 
# the commmon grid used between pileheight records. The group number along with Nx and Ny is written to a file. The second operation
# is to read the pileheight records agin and find the uniq_coords. Note that same sized grids would have same uniq_coords since the 
# downsampling aspect ratio is the same for them. The uniq_coord ( keys alone not actual co-ordinates) are written to a file.

if __name__ == '__main__':

    sample_Obj = UniqSampleGroups(2048)
    uniq_coords = [uniq_list for uniq_list, complete_list in itertools.groupby( sorted(sample_Obj.get_uniq_coord_list()) )]
# Write uniq_coords to uniq_coords.txt file

#    fp = open('uniq_coords_blah.txt','w')
    for coord_key in uniq_coords:
	 print coord_key
#        fp.write('%s\n' %coord_key)
#    fp.close()

# Group samples with same pileheight grid and write the data to file 
    fp = open('Nx_Ny_group_blah.txt','w')
    Map_coord = GridData(1)
    fp.write( str(Map_coord.X_start) + ' ' + str(Map_coord.X_end) + ' ' + str(Map_coord.Y_start) + ' ' + str(Map_coord.Y_end) + '\n' )
    uniq_grid = sample_Obj.get_uniq_Nx_Ny_dict()
    for key in uniq_grid:
        grid_list = uniq_grid[key]
        fp.write(str(key)+' ')
        fp.write(' '.join(str(grid_list).rstrip('\]').lstrip('\[').split(',')))
        fp.write('\n')
    fp.close()

