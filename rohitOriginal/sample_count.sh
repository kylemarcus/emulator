#! /bin/bash
nzsql -host netezza.ccr.buffalo.edu -u shivaswa -d rohit -pw Rohit123 -vON_ERROR_STOP=1 -t -A -F , <<-END
SELECT SC.sample,SS.count from sample_count as SC, (select * from (select sample,count(sample) over(partition by sample) from flowdata) as S group by sample,count) as SS where SC.sample=SS.sample order by SC.sample;
END
