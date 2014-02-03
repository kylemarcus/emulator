#! /bin/bash
nzsql -host netezza.ccr.buffalo.edu -u shivaswa -d rohit -pw Rohit123 -vON_ERROR_STOP=1 -t -A -F , <<-END
select count(*) from resamples;
select * from resamples order by resample;
END
