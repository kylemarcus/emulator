import re
import sys

inDir = "/panasas/scratch/kmarcus2/emulator/rohit/my_hadoop/"

resamples_dict = {}
resamples_weight_dict = {}
fp = open(inDir+"resamples.txt","r")
for line in fp:
  match = re.findall('([\d.]+)', line )
  resamples_dict[int(match[0])] =  match[1] + ' ' + match[2]  + ' ' + match[3] + ' ' + match[4]
  resamples_weight_dict[int(match[0])] = float(match[5])
fp.close()

#fp = sys.stdin
#fp = open("/panasas/scratch/kmarcus2/emulator/my_hadoop/resamples_used_2000/resamples_used.dat", 'r')
fp = open("/panasas/scratch/kmarcus2/emulator/my_hadoop/resamples_used_final.dat", 'r')
resamples = {}
for line in fp:
  resamples[int(line)] = 1;
fp.close()


phm = {}
fp = open(inDir+'montserrat_take2_vol_dir_bed_int.phm','r')
for line in fp:
  match = re.findall('[\d.]+',line)
  phm[int(match[0])] = match[1] + ' ' + match[2]
fp.close()

final = {}
sum = 0
#fp = open('final','r')
#fp = open("/panasas/scratch/kmarcus2/emulator/my_hadoop/reduce_output_2000/reduce_output.dat",'r')
fp = open("/panasas/scratch/kmarcus2/emulator/my_hadoop/reduce_output_final.dat",'r')
for line in fp:
  match = re.findall('[\d.]+',line)
  sum += resamples_weight_dict[int(match[0])]
  final[int(match[0])] = match[1] + ' ' + str(pow(float(match[1]),2))
fp.close()


fp = open('/scratch/hazard_map','w')
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
