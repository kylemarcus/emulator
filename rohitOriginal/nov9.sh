#! /bin/bash
nzsql -host netezza.ccr.buffalo.edu -u shivaswa -d rohit -pw Rohit123 -vON_ERROR_STOP=1 -t -A -F , <<-END
select P.X,P.Y,isnull(S1.count,0) from (select * from (select phm_id,count(sample) over(partition by phm_id) as count from (select phm_id,sample from phm_neighbours_nov11 group by phm_id,sample) as S) as SS group by phm_id,count) as S1 full join PHM as P on P.id=S1.phm_id order by P.id;
END
