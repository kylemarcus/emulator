#!/bin/bash
tablename=test_mean_and_variance$1
nzsql -host netezza -vON_ERROR_STOP=1 -t -A -F , <<-END
create table $tablename(SAMPLE int, UNIQ_ID int, RESAMPLE int, PHM_ID int, MEAN double precision, VARIANCE double precision) distribute on (sample,phm_id,resample);
END
