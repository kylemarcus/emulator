# ----------------------------------------------------------------------
#   MR-MPI = MapReduce-MPI library
#   http://www.cs.sandia.gov/~sjplimp/mapreduce.html
# -------------------------------------------------------------------------

# MapReduce emulator in Python
# Syntax: mpirun -np 8 python main.py /panasas/scratch/kmarcus2/emulator/mrmpi/down_sample_files

import sys
sys.path.append('.')
import subprocess
import reduce3
import reduce5
from mrmpi import mrmpi
try:
  import pypar
except:
  import pypar_serial as pypar


def mapEmulator(itask,filename,mr):
  
  # run the newEmulator with a downsample file
  cmd = 'python emulator.py ' + filename + ' | ./newEmulator'
  output = subprocess.check_output(cmd, shell=True)
  
  # get the list of output files
  o = output.strip().split('\n')
  for v in o:
    k = v.split(' ')[0]
    mr.add(k,v)
  
  print 'map done ' + str(itask)

def mapEmulatorOutputFiles(itask, filename, mr):

  # open emulator output file
  fp = open(filename, 'r')
  
  # add each line to the map
  for line in fp:
    key = line.strip().split(' ')[0]
    mr.add(key, line.strip())
  
  print 'map done ' + str(itask)

#def reduce3(itask, key,mvalue,mr):
def myreduce3(key,mvalue,mr):
  
  global reducePhase1Obj
  
  out = reducePhase1Obj.main(mvalue)
  
  if out != None:
    for line in out:
      key = '*' + line.split(' ')[0]
      mr.add(key,line)
    
  #print 'reduce3 done ' + str(key)

def myreduce5(key,mvalue,mr):
    
  global reducePhase2Obj
  
  if key[0] == '*':
    #fp = open('/panasas/scratch/kmarcus2/emulator/my_hadoop/reduce_output/'+str(key), 'w')
    
    out = reducePhase2Obj.main(mvalue)
    
    if out != None:
      fp = open('/scratch/reduce_output_'+str(key[1:]), 'w')
      for v in out:
        fp.write(v + '\n')
      fp.close()
    
    #print 'reduce5 done ' + str(key)







def reduceEmulator(key,mvalue,mr):
  
  cmd = 'cat '
  for f in mvalue:
    cmd += str(f) + ' '
  #cmd += '| python reduce3.py | python reduce4.py | python reduce5.py'
  cmd += '| python reduce3.py | python reduce5.py'
  output = subprocess.check_output(cmd, shell=True)
  o = output.split('\n')
  
  f = open('/panasas/scratch/kmarcus2/emulator/my_hadoop/reduce_output/'+str(key), 'w')
  for i in o:
    f.write(str(i))
  f.close()
  
  print 'reduce done ' + str(key)
  
  # reduce3.py
  # reduce4.py
  # reduce5.py


def mapPhmCount(itask,filename,mr):
  #print str(itask), filename
  l = map4.main(filename, 35)
  mr.add(itask+1, l)

def reducePhmCount(key,mvalue,mr):
  l = reduce4.main(mvalue[0])
  fn = '/panasas/scratch/kmarcus2/emulator/mrmpi/phm_count_files/phm_count_%d.txt' % key
  fp = open(fn,'w')
  for t in l:
    fp.write(str(t[0]) + ' ' + str(t[1]) + '\n')
  fp.close()
  

def mapTest(itask,filename,mr):
  l = [(itask,itask+1),(itask+2,itask+3)]
  mr.add(itask, l)

def reduceTest(key,mvalue,mr):
  print type(key)
  print mvalue
  print




# main program

nprocs = pypar.size()
me = pypar.rank()

if len(sys.argv) < 2:
  print "Syntax: emulator-mrmpi.py down_sample_files_dir"
  sys.exit()

files = sys.argv[1]

mr = mrmpi()
mr.verbosity(2)
mr.timer(1);

#mr2 = mr.copy(mr)

pypar.barrier()
tstart = pypar.time()

# each map function will get a different file name
#print "mr.map_file"
#rc = mr.map_file([files],0,0,0,mapEmulator)

print me, "mr.map_file"
rc = mr.map_file([files],0,0,0,mapEmulatorOutputFiles)

print me, "mr.collate after map"
rc = mr.collate()
print me, "collate rc: ", rc

#mr.kv_stats(2)
#mr.kmv_stats(2)

print me, "creating phase 1 reduce obj"
reducePhase1Obj = reduce3.ReducePhase1()

print me, "mr.reduce(myreduce3)"
rc = mr.reduce(myreduce3)
print me, "reduce rc: ", rc

print me, "mr.collate after reduce3"
rc = mr.collate()
print me, "collate rc: ", rc

#mr.kv_stats(2)
#mr.kmv_stats(2)

print me, "creating phase 2 reduce obj"
reducePhase2Obj = reduce5.ReducePhase2()

print me, "mr.reduce(myreduce5)"
rc = mr.reduce(myreduce5)
print me, "reduce rc: ", rc

print me, "mr.barrier"
pypar.barrier()
tstop = pypar.time()

#mr.sort_values(ncompare)
#count = [0,10,0]
#mr.map_mr(mr,output)

#mr.gather(1)
#mr.sort_values(ncompare)
#count = [0,10,1]
#mr.map_mr(mr,output)

mr.destroy()
#mr2.destroy()

# output

#if me == 0:
#  print "%d total words, %d unique words" % (nwords,nunique)
#  print "Time to process %d files on %d procs = %g (secs)" % \
#        (len(files),nprocs,tstop-tstart);
  
pypar.finalize()
