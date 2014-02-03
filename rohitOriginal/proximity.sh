#!/bin/bash
nzsql -host netezza.ccr.buffalo.edu -u shivaswa -d rohit -pw Rohit123 -vON_ERROR_STOP=1 -t -A -F , <<-END
select uniq_id-1,spatial_neighbour-1,uniq_id_neighbour-1,macro_neighbour-1 from proximity where sample=$1 order by spatialid;
END
