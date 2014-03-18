import sys
import pickle
#import itertools
from parser import *
from neighbor import *

def main():
    try:
	arg = int(sys.argv[1])
        return arg
    except ValueError:
	print '\nPlease enter the pileheight record number like 1,2,3, etc. To execute type "python downsample.py <integer number>\n'
	raise

if __name__ == '__main__':
    sample_number = main()

    sample_Obj = UniqSampleGroups(2048)
    my_dict = sample_Obj.get_sample_group_dict()
#    for key in my_dict:
#	print key, my_dict[key]

    group_index = sample_Obj.get_sample_group_index(sample_number)
    uniq_grid = sample_Obj.get_uniq_Nx_Ny_dict()
    print uniq_grid
    print 'sample group   ', sample_Obj.get_sample_group_index(sample_number)
    downsample_Obj = DownData(sample_number, sample_Obj.get_sample_group_index(sample_number) )
    downsampled_list = downsample_Obj.down_sample()

# Write downsampled data to the file. the first two columns are keys
    filename = 'down_sample_files/downsample%.6d.txt' %sample_number
    fp = open(filename ,'w')
    for sub_list in downsampled_list:
	for record in sub_list:
	    my_str = str(record)
	    new_str = (' '.join(re.split('\[|\'|\, |\]',my_str))).rstrip( ).lstrip( ) + '\n'
	    fp.write(new_str)

    fp.close()

