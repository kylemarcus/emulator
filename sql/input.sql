\echo 'start'

\set rohit /panasas/scratch/kmarcus2/emulator/input/rohit.txt
\echo :rohit
\set meta_data /panasas/scratch/kmarcus2/emulator/input/meta_data.txt
\echo :meta_data
\set macro_resamples /panasas/scratch/kmarcus2/emulator/input/macro_resamples.tmp
\echo :macro_resamples
\set phm /panasas/scratch/kmarcus2/emulator/input/montserrat_take2_vol_dir_bed_int.phm
\echo :phm

CREATE EXTERNAL TABLE NEWH(sample int,col int,row int,H double precision) using (dataobject(:rohit) Delimiter ',');
\echo 'creating external table NEWH from':rohit

CREATE TABLE TITAN AS SELECT * FROM newh;
\echo 'creating table TITAN from newh'

drop table newH;


CREATE EXTERNAL TABLE meta_ex(sample int, Nx int, Ny int, xstart double precision, xend double precision, ystart double precision, yend double precision,logvol double precision,direction double precision,basal double precision,internal double precision) using (dataobject(:meta_data) delimiter',');
\echo 'creaing external table meta_ex from':meta_data

CREATE TABLE meta as select * from meta_ex;
\echo 'creating table meta from meta_ex'

generate statistics on meta;

drop TABLE meta_ex;

create table temp as select titan.sample,titan.row,titan.col,titan.H,temp.X,temp.Y from titan inner join (select sample,row,col,((2*(row-1)+0.5)*(xend-xstart)/(2*Nx) + xstart) as X,((2*(col-1)+0.5)*(yend-ystart)/(2*Ny) + ystart) as Y from (select titan.sample,titan.row,titan.col,meta.xstart,meta.xend,meta.ystart,meta.yend,meta.Nx,meta.Ny from titan inner join meta on titan.sample=meta.sample) as joined) as temp on titan.sample=temp.sample and titan.row=temp.row and titan.col=temp.col order by sample;
\echo 'creating temp table'

DROP TABLE Titan;

ALTER TABLE temp rename to Titan;
\echo 'drop old Titan table and rename temp table to Titan'

GENERATE STATISTICS on TITAN;

create table downdata as select difference.sample, difference.dx, difference.dy,vol.volume,vol.geoflow_tiny from (select xdiff.sample,xdiff.dx,ydiff.dy from ((select part1.sample,(part1.X2-part2.X1) as dx from (select sample,titan.X as X2 from titan where row=2 and col=1) as part1 inner join (select sample,titan.X as X1 from titan where row=1 and col=1) as part2 on part1.sample=part2.sample)) as xdiff inner join ((select part3.sample,(part3.Y2-part4.Y1) as dy from (select sample,titan.Y as Y2 from titan where row=1 and col=2) as part3 inner join (select sample,titan.Y as Y1 from titan where row=1 and col=1) as part4 on part3.sample=part4.sample)) as ydiff on xdiff.sample=ydiff.sample) as difference inner join (select sample,pow(10,logvol) as volume, pow(pow(10,logvol),0.33333333333333)/1000 as GEOFLOW_TINY from meta) as vol on difference.sample=vol.sample order by sample;
\echo 'create table downdata'

/*
Nov 12,2012
The below line is commented out but can still be used. I have currently comenetd it and replaced it with another query right below it to capture more
points. The new query selects same number of points all over the map without distinction but certain samples have coarser grids and so the selection
of points is different from them. The final table of downsamle is a union of the two. This ensures sufficientnumber of points are selected within a
range of 100 metres even for coarser grids. The new query selects close to 6000 points per sample, but most can be skipped during emulator 
construction.

CREATE TABLE downsampled as select * from ( (select T.sample,T.row,T.col,T.X,T.Y,T.H from titan as T, downdata as D where D.sample=T.sample and T.row%2=0 and T.col%2=0 and T.H>3) union (select T.sample,T.row,T.col,T.X,T.Y,T.H from titan as T, downdata as D where D.sample=T.sample and T.row%2=0 and T.col%2=0 and T.H <= 3) ) as S where H > 0;
*/

CREATE TABLE downsampled as select * from ( (select T.sample,T.row,T.col,T.X,T.Y,T.H from titan as T, meta as M where M.sample=T.sample and T.row%3=0 and T.col%3=0 and M.Nx=320 and M.Ny=288) union (select T.sample,T.row,T.col,T.X,T.Y,T.H from titan as T, meta as M where M.sample=T.sample and T.row%2=0 and T.col%2=0 and M.Nx=160 and M.Ny=144) ) as S where H > 0;
\echo 'create table downsampled'

create table input_list as select sample, (logvol-(select min(logvol) from Meta))/( (select max(logvol) from meta) - (select min(logvol) from meta) ) as P1,(direction-(select min(direction) from Meta))/( (select max(direction) from meta) - (select min(direction) from meta) ) as P2,(basal-(select min(basal) from Meta))/( (select max(basal) from meta) - (select min(basal) from meta) ) as P3,(internal-(select min(internal) from Meta))/( (select max(internal) from meta) - (select min(internal) from meta) ) as P4 from meta;
\echo 'create table input_list'

create table copy_input_list as select * from input_list;
\echo 'create table copy_input_list'

CREATE TABLE macro_neighbour as select S.sample, S.neighbour,S.distance from (select P.sample as sample,C.sample as neighbour, sqrt( (P.p1-C.p1)*(P.p1-C.p1)+(P.p2-C.p2)*(P.p2-C.p2)+(P.p3-C.p3)*(P.p3-C.p3)+(P.p4-C.p4)*(P.p4-C.p4) ) as distance from input_list as P, copy_input_list as C) as S where S.distance<0.2 order by sample;
\echo 'create table macro_neighbour'

drop table copy_input_list;

CREATE TABLE flowdata as SELECT sample,rank() over(partition by sample order by row,col) as ID,X,Y,H from downsampled;
\echo 'create table flowdata'

DELETE from flowdata where H=0;

CREATE TABLE flowdata1 as select S2.uniq_id,S1.sample,S1.id,S1.x,S1.y,S1.h from (select rank() over(order by rank) as Uniq_id, rank from (select rank from (select *,rank() over(order by X,Y) from flowdata) as S group by rank) as SS) as S2,(select *,rank() over(order by X,Y) from flowdata) as S1 where S1.rank=S2.rank;
\echo 'create table flowdata1'

drop TABLE flowdata;

ALTER TABLE flowdata1 RENAME TO flowdata;

CREATE TABLE temporary as select M.sample,M.neighbour,F.ID,F.uniq_id,F.X,F.Y from macro_neighbour as M inner join flowdata as F on M.neighbour=F.sample; 
\echo 'create temporary table'

generate statistics on flowdata;

generate statistics on temporary;

/*
#<=100 denotes the distance along x and y axes. 200 is a preferred number fro the distance.
#rank<=150 indicates the number of neighbors you would like to store. It dictates the maximum size of the co-variance matrix.
*/

CREATE TABLE proximity as select SAMPLE,SPATIALID,MACRO_NEIGHBOUR,SPATIAL_NEIGHBOUR,UNIQ_ID,UNIQ_ID_NEIGHBOUR from (select *,rank() over(partition by sample,spatialid order by random) from (select F.sample,F.id as spatialid, T.neighbour as macro_neighbour, T.id as spatial_neighbour, F.uniq_id,T.uniq_id as uniq_id_neighbour,random() as random from flowdata as F, temporary as T where T.sample=F.sample and abs(F.X-T.X)<=100 and abs(F.Y-T.Y)<=100) as S) as SS where rank<=150;
\echo 'create table proximity'

CREATE TABLE proximity_details as select * from (select sample,spatialid,count(sample) over(partition by sample,spatialid) from proximity)  as S group by sample, spatialid, count;
\echo 'create table proximity_details'

drop table temporary;

CREATE EXTERNAL TABLE resamples_ex(resample int,P1 double precision,P2 double precision,P3 double precision,P4 double precision, w double precision) using (dataobject(:macro_resamples) delimiter ' ');
\echo 'create external table resamples_ex from':macro_resamples

CREATE TABLE resamples as select * from resamples_ex;
\echo 'create table resamples from resamples_ex'

drop TABLE resamples_ex;

CREATE TABLE maxmin_resamples as select max(p1) as maxp1, min(p1) as minp1, max(p2) as maxp2, min(p2) as minp2, max(p3) as maxp3, min(p3) as minp3, max(p4) as maxp4, min(p4) as minp4 from resamples;
\echo 'create table maxmin_resamples'

CREATE TABLE maxmin_samples as select max(logvol) as maxp1, min(logvol) as minp1, max(direction) as maxp2, min(direction) as minp2, max(basal) as maxp3, min(basal) as minp3, max(internal) as maxp4, min(internal) as minp4 from meta;
\echo 'create table maxmin_samples'

CREATE TABLE maxmin as select * from maxmin_resamples union select * from maxmin_samples;
\echo 'create table maxmin'

drop table maxmin_resamples,maxmin_samples;

CREATE TABLE minmax as select min(minp1) as minp1, max(maxp1) as maxp1, min(minp2) as minp2, max(maxp2) as maxp2, min(minp3) as minp3, max(maxp3) as maxp3, min(minp4) as minp4, max(maxp4) as maxp4 from maxmin;
\echo 'create table minmax'

CREATE TABLE samples_scaled as select sample,(S1.logvol-S2.minP1)/(maxP1 - minP1) as P1, (S1.direction-S2.minP2)/(maxP2 - minP2) as P2, (S1.basal-S2.minP3)/(maxP3 - minP3) as P3, (S1.internal-S2.minP4)/(maxP4 - minP4) as P4 from (select sample,logvol,direction,basal,internal from meta) as S1, (select * from minmax) as S2;
\echo 'create table samples_scaled'

drop TABLE maxmin;

CREATE TABLE resamples_scaled as select resample,(S1.p1-S2.minP1)/(maxP1 - minP1) as P1, (S1.p2-S2.minP2)/(maxP2 - minP2) as P2, (S1.p3-S2.minP3)/(maxP3 - minP3) as P3, (S1.p4-S2.minP4)/(maxP4 - minP4) as P4 from (select resample,P1,P2,P3,P4 from resamples) as S1, (select * from minmax) as S2;
\echo 'create table resamples_scaled'

CREATE TABLE Needed_resamples as select sample,resample, sqrt( (S1.P1-S2.P1)*(S1.P1-S2.P1) + (S1.P2-S2.P2)*(S1.P2-S2.P2) + (S1.P3-S2.P3)*(S1.P3-S2.P3) + (S1.P4-S2.P4)*(S1.P4-S2.P4) ) as distance from (select * from samples_scaled) as S1, (select * from resamples_scaled) as S2 where distance<.087;
\echo 'create table Needed_resamples'

CREATE TABLE new_scaled_resamples as select SS.rank as resample,S1.p1,S1.p2,S1.p3,S1.p4 from (select resample as old_resample, rank() over(order by resample) from (select resample from needed_resamples group by resample) as S) as SS, (select * from resamples_scaled) as S1 where SS.old_resample=S1.resample;
\echo 'create table new_scaled_resamples'

CREATE TABLE new_resamples as select SS.rank as resample,S1.p1,S1.p2,S1.p3,S1.p4,S1.w from (select resample as old_resample, rank() over(order by resample) from (select resample from needed_resamples group by resample) as S) as SS, (select * from resamples) as S1 where SS.rank=S1.resample;
\echo 'create talbe new_resamples'

drop TABLE resamples,resamples_scaled;

alter TABLE new_resamples RENAME TO resamples;

alter TABLE new_scaled_resamples RENAME TO resamples_scaled;

/* 0.12 dennotes the distance of search. You can change it to suit your need and your computational resource. */
CREATE TABLE resample_neighbours as select *,(1/distance)/sum(1/distance) over(partition by resample) as weight from (select resample,sample, sqrt( (S1.P1-S2.P1)*(S1.P1-S2.P1) + (S1.P2-S2.P2)*(S1.P2-S2.P2) + (S1.P3-S2.P3)*(S1.P3-S2.P3) + (S1.P4-S2.P4)*(S1.P4-S2.P4) ) as distance from (select * from samples_scaled) as S2, (select * from resamples_scaled) as S1 where distance<.12) as S3;
\echo 'create table resampled_neighbours'

drop table needed_resamples;

CREATE EXTERNAL TABLE phm_ex(id int,X double precision, Y double precision) using (dataobject(:phm) delimiter ' ');
\echo 'create external table phm_ex from':phm

CREATE TABLE phm as Select * from phm_ex;
\echo 'create table phm from phm_ex'

drop TABLE phm_ex;

CREATE TABLE uniq_coord as select uniq_id,X,Y from flowdata group by uniq_id,X,Y;
\echo 'create table uniq_coord'

CREATE TABLE temporary as select max(X) as maxx, min(x) as minx, max(y) as maxy, min(y) as miny from uniq_coord; 
\echo 'create table temporary'

CREATE TABLE phm_needed as SELECT P.id,P.x,P.y from phm as P,temporary as T where P.X <= T.maxx and P.X >= T.minx and P.Y <= T.maxy and P.Y >= T.miny;
\echo 'create table phm_needed'

drop TABLE tempoRARY;

CREATE TABLE phm_neighbours as select * from (select *,rank() over(partition by sample,phm_id order by distance) from (select S2.sample,S2.uniq_id,S2.phm_id,S4.spatialid,S2.distance from (select F.sample, F.id as spatialid, F.uniq_id, S1.phm_id,S1.distance from (select P.id as phm_id,U.uniq_id, sqrt((P.X-U.X)*(P.X-U.X)+(P.Y-U.Y)*(P.Y-U.Y)) as distance from uniq_coord as U, phm_needed as P where distance<=100 ) as S1 inner join flowdata as F on F.uniq_id=S1.uniq_id) as S2, (select S3.sample,S3.spatialid,F.uniq_id,S3.indicator from (select *,count/10 as indicator from proximity_details) as S3, flowdata as F where S3.sample=F.sample and S3.spatialid=F.id and S3.indicator>0) as S4 where S2.sample=S4.sample and S2.uniq_id=S4.uniq_id) as S) as SS where SS.rank<=3;
\echo 'create table phm_neighbours'

generate statistics on phm_neighbours;

/*
I have commnted the below 4 queries but they can be reused and must not be removed all together. While attempting to include more points through 
downsampling the below query which selects only those phm points that have 3 neighbours in spatial points was generating a jagged map with spikes at
reqular intervals. By not running the below queries the preliminary map which is map of the number of samples associated with phm points looks 
more smooth. So before running the emulator code on the cluster run the script "nov9.sh" and pipe out the output to a text file and view it on matlab.
"nov9.sh" will genetrate 3 columns, X co-ords, y co-ords and number of sample.

CREATE TABLE phm_neighbours_count as select *,(1/distance)/sum(1/distance) over(partition by sample,phm_id) as weight from (select S.sample,S.spatialid,S.uniq_id,S.phm_id,S.distance from (select *,count(phm_id) over(partition by sample,phm_id) from phm_neighbours) as S where count=3) as SS;

generate statistics on phm_neighbours_count;

drop table phm_neighbours;

alter TABLE phm_neighbours_count RENAME TO phm_neighbours;
*/

CREATE TABLE sample_count as select sample,sum(count) over(order by sample asc rows unbounded preceding exclude current row) from (select * from (select sample,count(sample) over(partition by sample) from proximity) as S group by sample,count) as SS;
\echo 'create table sample_count'

UPDATE sample_count set sum=0 where sum is null;

CREATE TABLE spatial_count as select *, sum(count) over(partition by sample order by spatialid asc rows unbounded preceding exclude current row) from (select SS.sample,SS.spatialid,F.H,SS.count from (select * from (select sample,spatialid,count(spatialid) over(partition by sample,spatialid) from proximity) as S group by S.spatialid,S.sample,S.count) as SS inner join flowdata as F on SS.sample=F.sample and SS.spatialid=F.id order by sample,spatialid) as SSS;
\echo 'create table spatial_count'

UPDATE spatial_count set sum=0 where sum is null;

CREATE TABLE weighted_phm(sample int, resample int, phm_id int, mean double precision, variance double precision) distribute on (phm_id,resample,sample);
\echo 'create table weighted_phm'

\echo 'done'
