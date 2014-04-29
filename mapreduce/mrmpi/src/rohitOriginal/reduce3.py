#!/usr/bin/env python
import sys
sys.path.append('.')
import re
from neighbor import *
from build_data_set import *
from parser import *

class Output:

  def __init__(self,key):
    self.key = key
    self.count = 0
    self.weight = 0
    self.mean = 0

  def get_count(self):
    return self.count

  def IncreaseCount(self):
    self.count += 1

  def UpdateMean(self,distance,mean,crit):
   if mean >= crit:
    self.mean += distance*mean
    self.weight += distance
    self.count += 1
   else:
    self.count += 1

  def get_weighted_mean(self):
    if self.weight == 0:
	return 0
    else:
	return self.mean/self.weight

  def get_weight(self):
    return self.weight

def ComputeDist(my_list1,my_list2):
  return pow ( (sum ( map ( (lambda x,y : (x-y)*(x-y)), my_list1, my_list2 ) )), 0.5 )


if __name__ == '__main__':
 Ndim = 4
 fp = open("/user/shivaswa/my_hadoop/resamples.txt","r")

 Resamples_dict = {}
 for i in range(Ndim):
    Resamples_dict[i+1] = []

 resamples_weight_dict = {}
 for line in fp:
   match = re.findall('([\d.\-e]+)', line )
   index = int(match[0])
   my_list = match[1:Ndim+1]
   for num,value in enumerate(my_list):
       Resamples_dict[num+1].append(float(value))
   resamples_weight_dict[int(match[0])] = float(match[Ndim+1])

#Scaling resamples data
 for key in Resamples_dict:
    my_list = Resamples_dict[key]
    min_list = min(my_list)
    max_list = max(my_list)
    Resamples_dict[key] = scale(my_list,max_list,min_list)


 Uncertain_List_Obj = Build_Sample_Object()
 Samples_dict = Uncertain_List_Obj.get_hashed_data()

 sample_neigh_count_dict = {}
 fp = open("/user/shivaswa/my_hadoop/resample_count.txt","r")
 for line in fp:
    match = re.search(r'(\d+)\s+(\d+)', line)
    sample_neigh_count_dict[int(match.group(1))] = int( match.group(2) )

 key = '0_0'
 count = 0
# ResPhm_Sam_dict = {}
 fp = sys.stdin
 for line in fp:
   try:
     match = re.findall( '([\d._\-e]+)', line)
   except:
     continue
   if len(match) != 3:
      continue
#   key = match[0] 
   new_key = match[0] 
   sample = int( match[1] )
   mean = float( match[2] )
   match = re.search(r'(\d+)_(\d+)',key)
   resample = int( match.group(1) )
   phm_id = int( match.group(2) )
   resample_list = []
   resample_list = [Resamples_dict[i+1][resample-1] for i in range(Ndim)]
   dist = ComputeDist(Samples_dict[sample],resample_list)

   if key != new_key:
      if count != 0:
        match = re.search(r'(\d+)_(\d+)',key)
        print match.group(2), '\t', match.group(1),' ', obj.get_weighted_mean()
      key = new_key
      obj = Output(key)
      obj.UpdateMean(dist,mean,0)
   else:
      count += 1
      obj.UpdateMean(dist,mean,0) 

 match = re.search(r'(\d+)_(\d+)',key)
 print match.group(2), '\t', match.group(1),' ', obj.get_weighted_mean()


# Below lines commented on Aug15, 2013
#   if key in ResPhm_Sam_dict:  
#      ResPhm_Sam_dict[key].UpdateMean(dist,mean,0)
#   else:
#      ResPhm_Sam_dict[key] = Output(key)
#      ResPhm_Sam_dict[key].UpdateMean(dist,mean,0)
  

# ***************       VERY old    ******************************************
#   if ( sample_neigh_count_dict[resample] == ResPhm_Sam_dict[key].get_count() ):
#     print phm_id, ' ', resample, ' ', ResPhm_Sam_dict[key].get_weighted_mean()
#     del ResPhm_Sam_dict[key]
# ****************    ********************  *************************************

# Below lines commented on Aug 15, 2013
# for key in ResPhm_Sam_dict:
#   match = re.search(r'(\d+)_(\d+)',key)
#   print match.group(2), '\t', match.group(1),' ', ResPhm_Sam_dict[key].get_weighted_mean()

#   del ResPhm_Sam_dict[key]

