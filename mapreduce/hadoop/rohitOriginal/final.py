import re
import sys


resamples_dict = {}
resamples_weight_dict = {}
fp = open("/user/shivaswa/my_hadoop/resamples.txt","r")
for line in fp:
   match = re.findall('([\d.]+)', line )
   resamples_dict[int(match[0])] =  match[1] + ' ' + match[2]  + ' ' + match[3] + ' ' + match[4]
   resamples_weight_dict[int(match[0])] = float(match[5])

fp = sys.stdin
resamples = {}
for line in fp:
  resamples[int(line)] = 1;
fp.close()


phm = {}
fp = open('montserrat_take2_vol_dir_bed_int.phm','r')
for line in fp:
  match = re.findall('[\d.]+',line)
  phm[int(match[0])] = match[1] + ' ' + match[2]
fp.close()

final = {}
sum = 0
fp = open('final','r')
for line in fp:
  match = re.findall('[\d.]+',line)
  sum += resamples_weight_dict[int(match[0])]
  final[int(match[0])] = match[1] + ' ' + str(pow(float(match[1]),2))
fp.close()


fp = open('hazard_map','w')
fp.write(str(len(resamples))+'\n')
fp.write(str(sum)+'\n')
for key in phm:
  try:
    fp.write(phm[key] + ' ' + final[key]+'\n')
  except:
    fp.write(phm[key] + ' ' +  '0' + ' ' + '0'+'\n')


#fp = open('resamples_used.txt','w')
for key in resamples:
  fp.write(resamples_dict[key]+'\n')
fp.close()
