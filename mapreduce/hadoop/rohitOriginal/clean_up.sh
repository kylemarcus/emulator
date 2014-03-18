#! /bin/bash

#cd /scratch
#rm -rf *
while read line
do
echo "node : " $line
ssh shivaswa@$line 'bash -s' < clear_scratch.sh
done
