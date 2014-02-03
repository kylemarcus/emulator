#! /bin/bash
echo "merging mean and variance...."
for i in `seq 1 $1`
do
echo $i
table1=mean_and_variance$i
table2=weighted_phm

nzsql -host netezza.ccr.buffalo.edu -vON_ERROR_STOP=1 -t -A -F , <<-END

create table mv_join as select M.sample,M.uniq_id,M.resample,P.phm_id,M.mean,M.variance,1/P.DISTANCE as distance from $table1 as M inner join phm_neighbours as P on  M.sample=P.sample and M.phm_id=P.phm_id and M.uniq_id=P.uniq_id;

create table mv_join1 as select sample,resample,phm_id,mean*distance/sum as mean, variance*distance/sum as variance from (select *,sum(distance) over(partition by sample,resample,phm_id) as sum from mv_join) as S distribute on (phm_id,sample,resample);

create table mv_join2 as select sample,resample,phm_id,sum(mean) over(partition by sample,resample,phm_id) as mean, sum(variance) over(partition by sample,resample,phm_id) as variance from mv_join1;

insert into $table2 select distinct * from mv_join2 ;

Generate statistics on $table2;

drop table mv_join,mv_join1,mv_join2;
END

done

nzsql -host netezza.ccr.buffalo.edu -vON_ERROR_STOP=1 -t -A -F , <<-END
CREATE TABLE mv_join_2(resample int, phm_id int, mean double precision) distribute on (phm_id,resample);
END

for i in `seq 1 5`
do

nzsql -host netezza.ccr.buffalo.edu -vON_ERROR_STOP=1 -t -A -F , <<-END

CREATE TABLE weighted_partition as  select * from weighted_phm where phm_id < $i*20000 and phm_id > ($i-1)*20000;

CREATE TABLE mv_join as select * from ( select W.sample,W.resample,W.phm_id,W.mean,W.variance,1/R.distance as distance from weighted_partition as W inner join resample_neighbours as R on W.sample=R.sample and W.resample=R.resample) as S;

create table mv_join_1 as select sample,resample,phm_id,mean*distance/sum as mean from (select *,sum(distance) over(partition by sample,resample,phm_id) as sum from mv_join) as S distribute on (phm_id,resample);

insert into mv_join_2 select distinct * from (select resample,phm_id,sum(mean) over(partition by resample,phm_id) as mean from mv_join_1) as S;

drop table weighted_partition,mv_join,mv_join_1;
END
done


nzsql -host netezza.ccr.buffalo.edu -vON_ERROR_STOP=1 -t -A -F , <<-END

CREATE TABLE final as select phm_id, sum(W) over(partition by phm_id) as W, sum(W2) over(partition by phm_id) as W2 from (select M.phm_id,round(floor(M.mean/0.2)/ceil((M.mean+0.0000001)/0.2))*R.w as W, round(floor(M.mean/0.2)/ceil((M.mean+0.0000001)/0.2))*R.w*R.w as W2 from mv_join_2 as M inner join resamples as R on M.resample=R.resample) as S;

CREATE TABLE final_view as select X,Y,isnull(w,0) as w,isnull(w2,0) as w2 from (select p.X,p.y,S.w,S.w2 from phm as p full join (select distinct * from final) as S on p.id=S.phm_id) as SS;

CREATE EXTERNAL TABLE resamples_used (P1 double precision, P2 double precision, P3 double precision, P4 double precision) using (dataobject('/home/shivaswa/final/resamples_used.txt') delimiter ' ');

insert into resamples_used select R.P1,R.p2,R.p3,R.p4 from (select resample from mv_join_2 group by resample) as S, resamples as R where S.resample=R.resample;

CREATE EXTERNAL TABLE view_file sameas final_view using (dataobject('/home/shivaswa/final/view_file.txt') delimiter ' ');

insert into view_file select * from final_view order by Y,X;

END

echo "NUMEBR of resample used...."

nzsql -host netezza.ccr.buffalo.edu -vON_ERROR_STOP=1 -t -A -F , <<-END
select count(*) from (select resample from mv_join_2 group by resample) as S;
END

echo "sum of weights of resamples used ....."
nzsql -host netezza.ccr.buffalo.edu -vON_ERROR_STOP=1 -t -A -F , <<-END
select sum(W) from (select R.w from (select resample from mv_join_2 group by resample) as S, resamples as R where S.resample=R.resample) as SS;
END
