# ----------------------------------------------------------------------
#   MR-MPI = MapReduce-MPI library
#   http://www.cs.sandia.gov/~sjplimp/mapreduce.html
# -------------------------------------------------------------------------

# MapReduce emulator in Python
# Syntax: mpirun -np 8 python main.py /panasas/scratch/kmarcus2/emulator/mrmpi/down_sample_files

import sys
import map4
import reduce4
from mrmpi import mrmpi
try:
  import pypar
except:
  import pypar_serial as pypar

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



# read a file
# for each word in file, emit key = word, value = NULL

def fileread(itask,mr):
  text = open(files[itask]).read()
  words = text.split()
  for word in words: mr.add(word,None)

# count word occurrence 
# emit key = word, value = # of multi-values

def sum(key,mvalue,mr):
  mr.add(key,len(mvalue))

# compare two counts
# order values by count, largest first

def ncompare(key1,key2):
  if key1 < key2: return 1
  elif key1 > key2: return -1
  else: return 0

# process a word and its count
# depending on flag, emit KV or print it, up to limit

def output(itask,key,value,mr):
  count[0] += 1
  if count[0] > count[1]: return
  if count[2]: print value,key
  else: mr.add(key,value)

# main program

nprocs = pypar.size()
me = pypar.rank()

if len(sys.argv) < 2:
  print "Syntax: phmCount-mrmpi.py down_sample_files_dir"
  sys.exit()
files = sys.argv[1]

mr = mrmpi()
mr.verbosity(2)
mr.timer(1);

pypar.barrier()
tstart = pypar.time()

# send the map function a list of files, in this case 1 dir
# each map function will get a different file name
#nwords = mr.map_file([files],0,0,0,mapEmulator)
nwords = mr.map_file([files],0,0,0,mapPhmCount)
mr.collate()
#nunique = mr.reduce(sum)
nunique = mr.reduce(reducePhmCount)

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

# output

#if me == 0:
#  print "%d total words, %d unique words" % (nwords,nunique)
#  print "Time to process %d files on %d procs = %g (secs)" % \
#        (len(files),nprocs,tstop-tstart);
  
pypar.finalize()
