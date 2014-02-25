#! /bin/bash
pipename=Netpipe$1
tablename=test_mean_and_variance$2
echo $pipename
nzload -host netezza.ccr.buffalo.edu -t $tablename -delim ' ' -df $pipename &
