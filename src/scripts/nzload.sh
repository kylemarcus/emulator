#! /bin/bash
pipename=Netpipe$1
tablename=mean_and_variance$2
echo $pipename
nzload -host netezza.ccr.buffalo.edu -t $tablename -delim ' ' -df $pipename &
