#!/usr/bin/env python
import sys
sys.path.append('.')
from reduce3 import Output
from neighbor import *
import re
import time

#def main(inputfp):
class ReducePhase2:
  
  def __init__(self):
    
    self.Ndim = 4
    fp = open("/panasas/scratch/shivaswa/emulator/my_hadoop/resamples.txt","r")
    self.resamples_weight_dict = {}
    for line in fp:
      match = re.findall('([\d.e-]+)', line )
      index = int(match[0])
      self.resamples_weight_dict[int(match[0])] = float(match[self.Ndim+1])
  #     print float(match[Ndim+1])

  #  obj = Output(1132)
    #fp = sys.stdin
  #  phm_res_dict = {}
  
  #@profile
  def main(self, inputfp):
    
    output = []
    resample_list = []
    phm_id = -1
    for line in inputfp:
      try:
        match = re.findall('([\d.\-e]+)',line )
        #match = line.split(' ')
      except:
        continue
      key = int(match[0])
  #    phm_id = int(match[0])
      resample = int(match[1])
      mean = float(match[2])
      resample_list.append(resample)
      weight = self.resamples_weight_dict[resample]

      if phm_id == key:
        obj.UpdateMean(weight,mean,0.2)

      else:
        if phm_id != -1:
          #print phm_id,' ',obj.get_weight()
          output += [str(phm_id) + ' ' + str(obj.get_weight())]
        phm_id = key
        obj = Output(phm_id)
        obj.UpdateMean(weight,mean,0.2)

    #print phm_id,' ',obj.get_weight()
    output += [str(phm_id) + ' ' + str(obj.get_weight())]
    
    #    if phm_id in phm_res_dict:
    #       phm_res_dict[phm_id].UpdateMean(mean,weight,0.2)
    #    else:
    #       phm_res_dict[phm_id] = Output(phm_id)
    #       phm_res_dict[phm_id].UpdateMean(mean,weight,0.2)


    #  for key in  phm_res_dict:
    #    print key,' ',phm_res_dict[key].get_weight()

    #filename = "/panasas/scratch/kmarcus2/emulator/my_hadoop/resamples_used/resample_" + str(phm_id)
    
    filename = "/scratch/resample_" + str(phm_id)
    fp = open(filename,"w")
    for resample in set(resample_list):
      fp.write(str(resample)+'\n')
    fp.close()
    
    return output

if __name__ == '__main__':
  #out = main(sys.stdin)
  #time1 = time.time()
  obj = ReducePhase2()
  #time2 = time.time()
  out = obj.main(sys.stdin)
  #time3 = time.time()
  #print time2 - time1
  #print time3 - time2
  for o in out:
    print o
