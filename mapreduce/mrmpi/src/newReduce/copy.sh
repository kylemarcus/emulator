#!/bin/sh

cd /scratch

cat reduce_output* > cat.reduce_output.$HOSTNAME.out &
#cat resample* > cat.resample.$HOSTNAME.out &
cat resample* | awk ' !x[$0]++' > cat.resample.$HOSTNAME.out &
wait

tar cfz reduce_output.$HOSTNAME.tar.gz cat.reduce_output.$HOSTNAME.out &
tar cfz resample.$HOSTNAME.tar.gz cat.resample.$HOSTNAME.out &
wait

cp reduce_output.$HOSTNAME.tar.gz /panasas/scratch/kmarcus2/emulator/my_hadoop/reduce_output_final &
cp resample.$HOSTNAME.tar.gz /panasas/scratch/kmarcus2/emulator/my_hadoop/resamples_used_final &
wait
