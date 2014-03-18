#!/usr/bin/env python
import sys
sys.path.append('.')
from reduce3 import Output
from neighbor import *
import re

if __name__ == '__main__':
  Ndim = 4

  fp = open("/panasas/scratch/shivaswa/emulator/my_hadoop/resamples.txt","r")
  resamples_weight_dict = {}
  for line in fp:
     match = re.findall('([\d.e-]+)', line )
     index = int(match[0])
     resamples_weight_dict[int(match[0])] = float(match[Ndim+1])
#     print float(match[Ndim+1])

#  obj = Output(1132)
  fp = sys.stdin
#  phm_res_dict = {}
  resample_list = []
  phm_id = -1
  for line in fp:
    try:
      match = re.findall('([\d.\-e]+)',line )
    except:
      continue
    key = int(match[0])
#    phm_id = int(match[0])
    resample = int(match[1])
    mean = float(match[2])
    resample_list.append(resample)
    weight = resamples_weight_dict[resample]

    if phm_id == key:
        obj.UpdateMean(weight,mean,0.2)

    else:
        if phm_id != -1:
          print phm_id,' ',obj.get_weight()
        phm_id = key
        obj = Output(phm_id)
        obj.UpdateMean(weight,mean,0.2)

  print phm_id,' ',obj.get_weight()
    
#    if phm_id in phm_res_dict:
#       phm_res_dict[phm_id].UpdateMean(mean,weight,0.2)
#    else:
#       phm_res_dict[phm_id] = Output(phm_id)
#       phm_res_dict[phm_id].UpdateMean(mean,weight,0.2)


#  for key in  phm_res_dict:
#    print key,' ',phm_res_dict[key].get_weight()

  filename = "/panasas/scratch/kmarcus2/emulator/my_hadoop/resamples_used/resample_" + str(phm_id)
  fp = open(filename,"w")
  for resample in set(resample_list):
    fp.write(str(resample)+'\n')
  fp.close()
