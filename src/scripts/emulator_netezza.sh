#!/bin/bash

# This is phase 1 in the emulator
#
# Before running this script, make sure you have the netezza module
# loaded and your environment vars set inside the module.
#
# Make sure you have the following files in some directory.
# Originally these files were in /home/shivaswa/final
# They have been moved to netezza:/nzscratch/emulator and /panasas/scratch/kmarcus2/emulator/input/
# ** You MUST specify these file locations in sql/input.sql
#
#    rohit.txt
#    meta_data.txt
#    montserrat_take2_vol_dir_bed_int.ph
#    macro_resamples.tmp
#    uncertain_input_list.txt <- this needs to be in extract dir, defined in cpp
#
# To run this script, you will have to change some of the variables. Change the following:
#  s - directory where sql files are located
#  o - directory where you want the ouput files to go
#  NZ_DATABASE - the database you want to use
# Also: for the first 2 nzsql command, you should specify your own database, USERNAME_db

NZ_DATABASE=testdb

s=/user/kmarcus2/emulator/sql
o=/panasas/scratch/kmarcus2/emulator/extract
sql=(${s}/input.sql ${s}/sample_count.sql ${s}/extract_phm.sql ${s}/extract_resamples.sql ${s}/uniq_coord.sql)
out=(${o}/input.out ${o}/sample_count.txt ${o}/phm.txt         ${o}/resamples.txt         ${o}/uniq_coords.txt)

nzsql -host netezza -d kmarcus2_db -c "DROP DATABASE $NZ_DATABASE"
nzsql -host netezza -d kmarcus2_db -c "CREATE DATABASE $NZ_DATABASE"

for i in ${!sql[*]}
do
	echo ${sql[$i]}
	nzsql -host netezza < ${sql[$i]} > ${out[$i]}
done
