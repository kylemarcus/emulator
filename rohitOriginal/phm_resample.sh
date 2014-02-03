#! /bin/bash
nzsql -host netezza.ccr.buffalo.edu -u shivaswa -d rohit -pw Rohit123 -vON_ERROR_STOP=1 -t -A -F , <<-END
select count(*) from (select resample from resample_neighbours where sample=$1) as S;
select resample-1 from resample_neighbours where sample=$1;
select isnull(count,0) from (select F.id,SS.count from (select * from (select spatialid,count(*) over(partition by sample,spatialid) from phm_neighbours where sample=$1) as S group by spatialid,count) as SS full join (select * from flowdata where sample=$1) as F on F.id=SS.spatialid) as SSS order by id;
select phm_id-1 from phm_neighbours where sample=$1 order by spatialid;
END




#select * from (select spatialid,count(*) over(partition by sample,spatialid) from phm_neighbours where sample=$1) as S group by spatialid,count order by spatialid;



