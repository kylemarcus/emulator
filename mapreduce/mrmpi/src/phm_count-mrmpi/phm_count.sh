#! /bin/bash

LIMIT=35
#for i in `seq 1 2048`
#do
filename=phm_count_$1.txt
cat uniq_phm.txt | python map4.py $1 $LIMIT | python reduce4.py > $filename
mv $filename /panasas/scratch/shivaswa/phm_count_files
#done

