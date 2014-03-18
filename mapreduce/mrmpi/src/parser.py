#!/usr/bin/env python

import sys
sys.path.append('.')
import re
import operator
import itertools

#generator function for evaluating modulo
def gen_mod(item_list,mod):
    for item in item_list:
        if item % mod == 0:
            yield item

def hash_key(sample,row,col): return ( str(sample) + '_' +  str(row) + '_' + str(col) )
def uniq_coord_key(group_index,row,col) : return ( str(group_index) + '_' +  str(row) + '_' + str(col) )

def scale(my_list, max_list = 0, min_list = 0):
    if max_list == min_list:
        max_list = float(max(my_list))
        min_list = float(min(my_list))
        print max_list, min_list
    return [ float((item-min_list)/(max_list-min_list)) for item in my_list]
    

def transpose(data_list):
    my_list = [list(tup) for tup in zip(*data_list)]
    return my_list


# downsampl_aspect_ratio allows you to control the extent of downsampling. Evry grid cell whose modulo with this number returns 
# zero is admitted in the downsample list, rest of the points are discarded. Bigger the aspect ratio more the points discarded.
def downsample_aspect_ratio(num): 
    if num > 384:
        aspect_ratio = 5
    elif num <= 384 and num  > 288:
        aspect_ratio = 4
    elif num <= 288 and num  >= 192:
        aspect_ratio = 3
    else:
        aspect_ratio = 2
    print "  aspect ratio = " , aspect_ratio
    return aspect_ratio


class GridData:
    def __init__(self,sample):

        self.filename = '/panasas/scratch/shivaswa/pileheight/pileheightrecord.%.6d' % sample
        fp = open(self.filename,'r')
        line1 = re.search(r'(\d+).+\s+([\d.]+),\s+([\d.]+)', fp.readline() )
        line2 = re.search(r'(\d+).+\s+([\d.]+),\s+([\d.]+)', fp.readline() )

        self.Nx = int(line1.group(1))
        self.Ny = int(line2.group(1))
        self.X_start = float(line1.group(2))
        self.X_end = float(line1.group(3))
        self.Y_start = float(line2.group(2))
        self.Y_end = float(line2.group(3))
        self.sample = sample
        fp.close()

    def get_filename(self): return self.filename
    def get_Nx(self): return self.Nx
    def get_Ny(self): return self.Ny
    def get_sample(self): return self.sample

# Nx is the number of columns in pileightrecord files
# Ny is the number of rows in pileightrecord files
class UniqSampleGroups:
    def __init__(self,num_of_samples, grid_size_list=[]):
        self.Nx_Ny_list = grid_size_list

        if len(grid_size_list) != num_of_samples:
            self.Nx_Ny_list = [ self.__Read_Nx_Ny(i+1) for i in range(num_of_samples) ]

        self.sample_group = {}
        self.__create_sample_groups()

    def __Read_Nx_Ny(self,num):
        obj = GridData(num)
        return [num,obj.get_Nx(),obj.get_Ny()]

    def __gen_group(self,my_list):
        for sub_list in my_list:
            sample_group = []
            for sub in sub_list:
                sample_group.append(sub[0])
            yield sample_group

    def __increment(self,num):
        num[0] += 1
        return num[0]

    def __create_sample_groups(self):
        sorted_list = sorted ( self.Nx_Ny_list, key = operator.itemgetter(1,2) )
        grouped_list = [list(group_list) for key_list, group_list  in itertools.groupby( sorted_list, key = operator.itemgetter(1,2) )]
        count = [0]
        self.sample_group = dict( (self.__increment(count), list(iter))  for iter in self.__gen_group(grouped_list) )
        #OR sample_group = { (increment(count) : list(iter))  for iter in gen_group(grouped_list) }

    def get_sample_group_dict(self):
        return self.sample_group

    def get_uniq_Nx_Ny_dict(self):
#       sample_group = create_sample_groups()
        uniq_Nx_Ny = {}
        for key in self.sample_group:
            Nx = self.Nx_Ny_list[ ( self.sample_group[key] )[0] - 1 ][1]
            Ny = self.Nx_Ny_list[ ( self.sample_group[key] )[0] - 1 ][2]
            uniq_Nx_Ny[key] = [Nx,Ny]
        return uniq_Nx_Ny

    def get_sample_group_index(self,sample):
        try:
            1/len(self.sample_group)
        except ZeroDivisionError:
            print ' Invoke function "__create_sample_groups()" before calling this funtion '

        for key in self.sample_group:
            try:
                group_index = self.sample_group[key].index(sample)
                break
            except:
                print 'BLAH'
                continue
        return key

    def get_uniq_coord_list(self):
        uniq_Nx_Ny  = self.get_uniq_Nx_Ny_dict()
        uniq_coord_list = []
        for key in uniq_Nx_Ny:
            my_list = []
            Nx = uniq_Nx_Ny[key][0]
            Ny = uniq_Nx_Ny[key][1]
            mod_row = downsample_aspect_ratio(Ny)
            mod_col = downsample_aspect_ratio(Nx)
            row_list = [i+1 for i in range(Ny)]
            col_list = [i+1 for i in range(Nx)]
            my_list = [str(key)+'_'+str(row)+'_'+str(col) for col in gen_mod(col_list,mod_col) for row in gen_mod(row_list,mod_row)]
            uniq_coord_list = uniq_coord_list + my_list
        return uniq_coord_list

class DownData:

    def __init__(self,sample,group_index):
        self.grid_info = GridData(sample)
        self.group_index = group_index

    def down_sample(self):
        sample = self.grid_info.get_sample()
        Nx = self.grid_info.get_Nx()    
        Ny = self.grid_info.get_Ny()    
        mod_row = downsample_aspect_ratio(Ny)
        mod_col = downsample_aspect_ratio(Nx)
        print " Nx = " , Nx , "    Ny = " , Ny
        print "mod_row = ", mod_row, "   mod_col = ",mod_col
        row_list = [i+1 for i in range(Ny)]
        col_list = [i+1 for i in range(Nx)]

        fp = open(self.grid_info.get_filename(),'r')
        [fp.readline() for i in range(3)]     #Skipping first 3 lines of pileheight record which contains metadata
        down_list = [list(self.__sift_records(list(pile_list)[0],list(pile_list)[1],mod_col)) for pile_list in self.__gen_line(fp,row_list,mod_row)]
        fp.close()
        return down_list

    def __gen_line(self,fp,row_list,mod):
        for row in row_list:
            read_line = fp.readline()   
            if row % mod == 0:
                yield [row,list(read_line.rstrip('\n').split(' '))]

    def __sift_records(self,row,pile_list,mod):
        sample = self.grid_info.get_sample()
        count = 0
        for pile in pile_list:
            float_pile = float(pile)
            count += 1
            if count % mod == 0 and float_pile > 0:
                yield [sample,uniq_coord_key(self.group_index,row,count), float_pile]
#               yield [hash_key(sample,row,count), uniq_coord_key(self.group_index,row,count), float_pile]


# This reads the data from a file and stores it as a list of list. It can also store unstructured dataset. 
# fp stands for file pointer and should be passed when pointing to the data.
class DataRead():

    def __init__(self):
        self.num_row = 0
        self.num_col = 0

    def load_from_file(self,fp):
        my_list = []
        for i,line in enumerate(fp):
            match = re.findall('[\d\.\_\-e]+',line)
            my_list.append([num for num in self.__gen_number(match)])
            count = i
        self.num_row = count
        self.num_col = len(my_list[0])
        return my_list

    def rows(self): return self.num_row

    def cols(self): return self.num_col

    def __gen_number(self,match):
        for str in match:
            if re.search(r'\_',str):
                yield str
            elif re.search(r'[\.\-e]',str):
                yield float(str)
            else:
                yield int(str)


def get_max_min(fp):
    read_data = DataRead()
    data = read_data.load_from_file(fp)
    tran_data = transpose(data)
    my_dict = dict( (key, [max(sub_list), min(sub_list)]) for key,sub_list in enumerate(tran_data) )
    return my_dict

def gen_scaled_data(my_list,max_min_dict):
    for key,sub_list in enumerate(my_list):
        if key == 0: yield sub_list
        else:
            scaled_list =  scale(sub_list, max_min_dict[key][0], max_min_dict[key][1])
            yield scaled_list

def print_key_value(key,my_list):
    for item in my_list:
        print item,key

def gen_key_value(fp):
    for line in fp:
        match = re.findall('[\d._]+',line)
        yield [int(match[0]), int(match[1])]

def gen_reverse(my_list):
    for item in my_list:
        yield item[1]
