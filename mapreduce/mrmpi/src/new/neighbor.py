#!/usr/bin/env python

import sys
sys.path.append('.')
from parser import *

def gen_zip_list(my_list):
    for zipped_tuple in zip(*my_list):
        yield list(zipped_tuple)


def disp_dictionary(my_dict):
    for key in my_dict:
        print 'key = %s' % key , '  data = ' , my_dict[key]


class StructuredDataSet:

    def __init__(self,two_dim_list=[]):
        self.data_list = two_dim_list
        self.data_size = len(two_dim_list[0])
        self.num_param = len(two_dim_list)
        self.hash_dict = dict
        self.hash_column = -1

# Scaling the data set replaces original data set with the scaled data set and does not maintain a copy.
    def scale(self,hash_column = -1):
#        scaled_data = [scale(column, max(column) ,min(column)) for column in self.data_list]
        scaled_data = [scale(column, max(column), min(column)) if hash_column != num else column for num, column in enumerate(self.data_list)]
        self.data_list = scaled_data

# The default hash key (if no column is identified as key) starts from 1 and not 0.
    def hash_map(self, hash_column = -1):
        self.hash_column = hash_column
        param_key_list = map( lambda x : x + 1, range(self.num_param) ) if hash_column == -1 else self.data_list[hash_column]
        data_set_without_key = [sub_list for i,sub_list in enumerate(self.data_list) if i != hash_column]
        key_list = map( (lambda x:x+1),range(self.data_size) ) if hash_column == -1 else self.data_list[hash_column]
        self.hash_dict = dict( (my_tuple[0], my_tuple[1]) for my_tuple in zip( key_list, list( gen_zip_list(data_set_without_key) ) ) )

#    def disp_key_value(self):
#        for key in self.hash_dict:
#            print 'key = %s' % key , '  data = ' , self.hash_dict[key]

    def get_hashed_data(self): return self.hash_dict

    def get_key_list(self) : return set([key for key in self.hash_dict])

    def get_parameter(self,key): return self.hash_dict[key]

#    def find_neighbors(self,key_value_dict,dist):
#       neigh_dict = dict( ( key, self.__gen_neigh_key_list(self.hash_dict[key] ,key_value_dict,dist,key) ) for key in self.hash_dict)
#       return neigh_dict

    def find_neighbors(self,key_value_dict,dist):
        neigh_dict = dict( ( key, self.__gen_neigh_key_list(self.hash_dict[key] ,key_value_dict,dist) ) for key in self.hash_dict)
        return neigh_dict

    def __gen_neigh_key_list(self, current_list, neigh_data, dist):
        key_list = [key for key in neigh_data if pow(sum(map(( lambda x,y : (x-y)*(x-y)), current_list, neigh_data[key])), 0.5) <= dist]
        return key_list

#    def __gen_neigh_key_list(self, current_list, neigh_data, dist, key):
#       key_list = [key for key in neigh_data if pow(sum(map(( lambda x,y : (x-y)*(x-y)), current_list, neigh_data[key])), 0.5) <= dist]
#       return key_list

    def get_neighbors(self,param_list,dist):
        return self.__gen_neigh_key_list(param_list,self.hash_dict,dist)    

    def gen_find_neighbors_with_dist(self,key_value_dict,dist):
        for neigh_key in key_value_dict:
            for my_key in self.hash_dict:
                my_dist = pow(sum(map(( lambda x,y : (x-y)*(x-y)), key_value_dict[neigh_key], self.hash_dict[my_key])), 0.5)
                if my_dist <= dist:
                    yield [neigh_key, my_key, my_dist]


    def key_intersect(self,key_set): return  self.get_key_set().intersection(set(key_set))
                
