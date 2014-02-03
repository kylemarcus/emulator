#!/bin/bash
tablename=TITAN
pipename=mypipe
nzload -host netezza.ccr.buffalo.edu -u shivaswa -pw Rohit123 -db rohit -t $tablename -delim ',' -df $pipename &
