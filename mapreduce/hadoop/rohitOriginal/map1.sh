#! /bin/bash

cd ~/my_hadoop
filename="POP_"`hostname`
date > $filename
python map1.py $filename
