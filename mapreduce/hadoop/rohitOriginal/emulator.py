#!/usr/bin/env python
import sys
sys.path.append('.')
from parser import *
from build_data_set import *
import time
import numpy as np
import numpy.random as rnd


def get_int_list(my_list):
    return map( (lambda x : int(x)), my_list )

def get_list_intersection(my_list1,my_list2):
    return list(set(my_list1).intersection(set(my_list2)))

def extract_downsampled_data(sample):
    filename = '/user/shivaswa/my_hadoop/down_sample_files/downsample%0.6d.txt' %sample
    fp = open(filename,'r')
    down_data_file = DataRead()
    down_data = down_data_file.load_from_file(fp)
    downsample = StructuredDataSet( transpose(down_data)[1:3] )
    downsample.hash_map(0)
    fp.close()
    return downsample
#    return downsample.get_hashed_data()

def gen_coords_pileheight(common_key_list,xy_coord_dict, downsample_dict):
    for key in common_key_list:
	yield xy_coord_dict[key]+downsample_dict[key]

def gen_uniq_intersection_key(my_list1,my_list2):
    for common_uniq_key in get_list_intersection(my_list1,my_list2):
	yield common_uniq_key	    

def gen_randomX(my_size,two_d_list,cut):
#    for num in (list(np.ceil(my_size*(rnd.random(50))) ))[0:cut]:
    for num in list(set(np.ceil( (my_size*(rnd.random(500))).tolist() ) ))[0:cut]:
        yield two_d_list[int(num)]

if __name__ == '__main__':

    MACRO_NEIGH_LIMIT = 320
    PHM_NEIGH_LIMIT = 35

#Read uncertain input list
    fp = open('/user/shivaswa/my_hadoop/uncertain_input_list.txt','r')
    skip_lines = int(re.search(r'\d+',fp.readline()).group())

    junk = [fp.readline() for i in range(skip_lines+2)]

    unc_inp_file = DataRead()
    sample_data = unc_inp_file.load_from_file(fp)
    tran_sample_data = transpose(sample_data)

    for col in tran_sample_data:
        print max(col), min(col) 
    input_parameter = StructuredDataSet(tran_sample_data)
    input_parameter.hash_map()
    input_parameter_dict = input_parameter.get_hashed_data()

#Set up sample_sample_neigh_dict
    sample_obj = Build_Sample_Object()
    sample_key_dict = sample_obj.get_hashed_data()
    sample_sample_neigh_dict = sample_obj.find_neighbors(sample_key_dict,0.2)
#    disp_dictionary(sample_sample_neigh_dict)
#NoTE : uiq_uniq.txt doesn't yet exist. To create it you need to run map3.py on hadoop and then merge all the output files. 
# Right now map3.py prints out the output to stdout. I'm yet to decide whether to let haddop write t tpo files lke the reduce step 
# does or my self write to files.In either case it should be noted that it has no reduce step. The merged file will be called
# uniq_uniq.txt..

#Set up uniq_uniq_dict
    uniq_uniq_neigh_dict = {}
    fp = open('/user/shivaswa/my_hadoop/uniq_uniq.txt','r')
    for line in fp:
        match = re.findall('([\d_]+)',line)
        uniq_uniq_neigh_dict[match.pop(0)] = match
    fp.close()

# Set up sample_resample_dict
    sample_resample_neigh_dict = {}
    fp = open('/user/shivaswa/my_hadoop/sample_resample.txt','r')
    for line in fp:
        match = re.findall('([\d]+)',line)
        sample_resample_neigh_dict[int(match.pop(0))] = get_int_list(match)
    fp.close()

# Set up uniq_phm_neigh_dict
    uniq_phm_neigh_dict = {}
    fp = open('/user/shivaswa/my_hadoop/uniq_phm.txt','r')
    for line in fp:
        match = re.findall('([\d_]+)',line)
	key = match.pop(0)
        uniq_phm_neigh_dict[key] = get_int_list(match)
    fp.close()
#    disp_dictionary(uniq_phm_neigh_dict)

    uniq_coords_obj = Build_Uniq_Coords_Object()
    xy_coord_dict =  uniq_coords_obj.get_hashed_data()

    sample = 0
#07/25/13 - Commenting the below line to change map function from hadoop based to mpi scheduler based
#    fp = sys.stdin
    arg_sample = int(sys.argv[1])
    filename = 'down_sample_files/downsample%06d.txt'%arg_sample
    fp = open(filename,'r')

    for line in fp:
	match = re.findall('([\d_.\-e]+)',line)    		
#	my_sample_num = int(match(0))
	my_sample_num = int(match[0])
	my_uniq_coords_key = match[1]

	my_pileheight = float( match[2] )
    	if sample != my_sample_num:
	    sample = my_sample_num
	    print 'SAMPLE'
	    print sample
	    print 'STOP'
	    my_sample_neigh = [neigh for neigh in sample_sample_neigh_dict[sample] if neigh is not sample]
	    my_sample_neigh = [sample] + my_sample_neigh
	    print 'RESAMPLES'

	    try:
	       sample_resample_neigh_list = sample_resample_neigh_dict[sample]
	    except KeyError:
	       sys.exit()
	    MY_LIMIT = 320 if len(sample_resample_neigh_list) > MACRO_NEIGH_LIMIT else sample_resample_neigh_list

	    my_resample_neigh = [neigh for num,neigh in enumerate(sample_resample_neigh_list) if num < MY_LIMIT ]
	    print ' '.join( (re.split(',' , str(my_resample_neigh).rstrip('\]').lstrip('\[')) ) )
	    print 'STOP'
	    sample_downsample_dict = dict( (sam_neigh,extract_downsampled_data(sam_neigh)) for sam_neigh in my_sample_neigh )
	    # sample_downsample_dict is a dictionary of StructuredDataSet objects. Each object has downsample information extracted
	    # from downsmapled file.

	my_uniq_coord_neigh = uniq_uniq_neigh_dict[my_uniq_coords_key]
	X = []

	try:
          my_X = xy_coord_dict[my_uniq_coords_key] + input_parameter_dict[sample] + [my_pileheight]
	except:
	  print 'my_uniq_coords_key does not exist. You might want to check files -- uniq_uniq.txt or uniq_phm.txt.'

	for item in my_sample_neigh:
	    param_set = input_parameter_dict[item]
            downsample_obj = sample_downsample_dict[item]
	    iter_common_key = gen_uniq_intersection_key( my_uniq_coord_neigh, downsample_obj.get_key_list() )
	    d = [xy_coord_pile[0:2]+param_set+[xy_coord_pile[2]] for xy_coord_pile in gen_coords_pileheight(iter_common_key,xy_coord_dict,downsample_obj.get_hashed_data())]
	    X = X+d
	    for index,X_line in enumerate(X):
		if X_line == my_X:
		    X.pop(index)
        if ( len(X) > 200 ):
	    print 'NEXT'
	    print my_uniq_coords_key
	    print ' '.join(re.split('\'|, ', str(my_X).rstrip('\]').lstrip('\[')))
	    for my_list in gen_randomX(len(X)-1,X,199):
	        print ' '.join(re.split('\'|, ', str(my_list).rstrip('\]').lstrip('\[')))
	else:
	    if len(X) == 0:
		print 'STTTTTTTTTTTTTTTTTTTTTOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP'
	    print 'NEXT'
	    print my_uniq_coords_key
	    print ' '.join(re.split('\'|, ', str(my_X).rstrip('\]').lstrip('\[')))
	    for my_list in X:
		print ' '.join(re.split('\'|, ', str(my_list).rstrip('\]').lstrip('\[')))
	print 'STOP'

	print 'PHM'
	phm_neighbors_list = uniq_phm_neigh_dict[my_uniq_coords_key];
	phm_limit = PHM_NEIGH_LIMIT if len(phm_neighbors_list) > PHM_NEIGH_LIMIT else len(phm_neighbors_list)
        print ' '.join(re.split('\, ', str( phm_neighbors_list[0:phm_limit] ).rstrip('\]').lstrip('\[') ) )
#        print ' '.join(re.split('\, ', str( (uniq_phm_neigh_dict[my_uniq_coords_key])[0:35] ).rstrip('\]').lstrip('\[') ) )
	print 'STOP'
#	print ' '.join(re.split('\[|\]|, ' ,str(uniq_phm_neigh_dict[my_uniq_coords_key]) ) )


    fp.close()
