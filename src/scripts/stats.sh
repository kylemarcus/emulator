#!/bin/bash
tablename=mean_and_variance$1
echo GENERATING STATISTICS - Message from stats.sh 
echo "TABLE" $1
nzsql -host netezza.ccr.buffalo.edu -vON_ERROR_STOP=1 -t -A -F , <<-END
GENERATE STATISTICS ON $tablename;
END
