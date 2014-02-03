#! /bin/bash
nzsql -host netezza.ccr.buffalo.edu -vON_ERROR_STOP=1 -t -A -F , <<-END
select H,COUNT from spatial_count where sample <= 2048 order by sample,spatialid;
END
