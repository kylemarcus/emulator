CREATE TABLE titan(sample int, col int, row int, H double PRECISION);
LOAD DATA LOCAL INFILE '/keith1/data/users/kmarcus2/sql-input/phData.txt' 
  INTO TABLE titan
  FIELDS TERMINATED BY ',';

CREATE TABLE meta(sample int, Nx int, Ny int, xstart double PRECISION, xend double PRECISION, ystart double PRECISION, yend double PRECISION,logvol double PRECISION,direction double PRECISION,basal double PRECISION,internal double PRECISION);
LOAD DATA LOCAL INFILE '/keith1/data/users/kmarcus2/sql-input/phMetaData.txt' 
  INTO TABLE meta
  FIELDS TERMINATED BY ',';

CREATE TABLE TEMP AS
SELECT titan.sample,
       titan.ROW,
       titan.col,
       titan.H,
       TEMP.X,
       TEMP.Y
FROM titan
INNER JOIN
  (SELECT sample,
          ROW,
          col,
          ((2*(ROW-1)+0.5)*(xend-xstart)/(2*Nx) + xstart) AS X,
          ((2*(col-1)+0.5)*(yend-ystart)/(2*Ny) + ystart) AS Y
   FROM
     (SELECT titan.sample,
             titan.ROW,
             titan.col,
             meta.xstart,
             meta.xend,
             meta.ystart,
             meta.yend,
             meta.Nx,
             meta.Ny
      FROM titan
      INNER JOIN meta ON titan.sample=meta.sample) AS joined) AS TEMP ON titan.sample=TEMP.sample
AND titan.ROW=TEMP.ROW
AND titan.col=TEMP.col
ORDER BY sample;

DROP TABLE titan;

ALTER TABLE TEMP RENAME TO Titan;

CREATE TABLE downdata AS
SELECT difference.sample,
       difference.dx,
       difference.dy,
       vol.volume,
       vol.geoflow_tiny
FROM
  (SELECT xdiff.sample,
          xdiff.dx,
          ydiff.dy
   FROM (
           (SELECT part1.sample,
                   (part1.X2-part2.X1) AS dx
            FROM
              (SELECT sample,
                      titan.X AS X2
               FROM titan
               WHERE ROW=2
                 AND col=1) AS part1
            INNER JOIN
              (SELECT sample,
                      titan.X AS X1
               FROM titan
               WHERE ROW=1
                 AND col=1) AS part2 ON part1.sample=part2.sample)) AS xdiff
   INNER JOIN (
                 (SELECT part3.sample,
                         (part3.Y2-part4.Y1) AS dy
                  FROM
                    (SELECT sample,
                            titan.Y AS Y2
                     FROM titan
                     WHERE ROW=1
                       AND col=2) AS part3
                  INNER JOIN
                    (SELECT sample,
                            titan.Y AS Y1
                     FROM titan
                     WHERE ROW=1
                       AND col=1) AS part4 ON part3.sample=part4.sample)) AS ydiff ON xdiff.sample=ydiff.sample) AS difference
INNER JOIN
  (SELECT sample,
          pow(10,logvol) AS volume,
          pow(pow(10,logvol),0.33333333333333)/1000 AS GEOFLOW_TINY
   FROM meta) AS vol ON difference.sample=vol.sample
ORDER BY sample;

/*
Nov 12,2012
The below line is commented out but can still be used. I have currently comenetd it and replaced it with another query right below it to capture more
points. The new query selects same number of points all over the map without distinction but certain samples have coarser grids and so the selection
of points is different from them. The final table of downsamle is a union of the two. This ensures sufficientnumber of points are selected within a
range of 100 metres even for coarser grids. The new query selects close to 6000 points per sample, but most can be skipped during emulator
construction.

CREATE TABLE downsampled as select * from ( (select T.sample,T.row,T.col,T.X,T.Y,T.H from titan as T, downdata as D where D.sample=T.sample and T.row%2=0 and T.col%2=0 and T.H>3) union (select T.sample,T.row,T.col,T.X,T.Y,T.H from titan as T, downdata as D where D.sample=T.sample and T.row%2=0 and T.col%2=0 and T.H <= 3) ) as S where H > 0;
*/

CREATE TABLE downsampled AS
SELECT *
FROM (
        (SELECT T.sample,
                T.ROW,
                T.col,
                T.X,
                T.Y,
                T.H
         FROM titan AS T,
              meta AS M
         WHERE M.sample=T.sample
           AND T.ROW%3=0
           AND T.col%3=0
           AND M.Nx=320
           AND M.Ny=288)
      UNION
        (SELECT T.sample,
                T.ROW,
                T.col,
                T.X,
                T.Y,
                T.H
         FROM titan AS T,
              meta AS M
         WHERE M.sample=T.sample
           AND T.ROW%2=0
           AND T.col%2=0
           AND M.Nx=160
           AND M.Ny=144)) AS S
WHERE H > 0;

CREATE TABLE input_list AS
SELECT sample,
       (logvol-
          (SELECT min(logvol)
           FROM Meta))/(
                          (SELECT max(logvol)
                           FROM meta) -
                          (SELECT min(logvol)
                           FROM meta)) AS P1,
       (direction-
          (SELECT min(direction)
           FROM Meta))/(
                          (SELECT max(direction)
                           FROM meta) -
                          (SELECT min(direction)
                           FROM meta)) AS P2,
       (basal-
          (SELECT min(basal)
           FROM Meta))/(
                          (SELECT max(basal)
                           FROM meta) -
                          (SELECT min(basal)
                           FROM meta)) AS P3,
       (internal-
          (SELECT min(internal)
           FROM Meta))/(
                          (SELECT max(internal)
                           FROM meta) -
                          (SELECT min(internal)
                           FROM meta)) AS P4
FROM meta;

CREATE TABLE copy_input_list AS
SELECT *
FROM input_list;

CREATE TABLE macro_neighbour AS
SELECT S.sample,
       S.neighbour,
       S.distance
FROM
  (SELECT P.sample AS sample,
          C.sample AS neighbour,
          sqrt((P.p1-C.p1)*(P.p1-C.p1)+(P.p2-C.p2)*(P.p2-C.p2)+(P.p3-C.p3)*(P.p3-C.p3)+(P.p4-C.p4)*(P.p4-C.p4)) AS distance
   FROM input_list AS P,
        copy_input_list AS C) AS S
WHERE S.distance<0.2
ORDER BY sample;

DROP TABLE copy_input_list;

CREATE TABLE flowdata AS
SELECT sample,
       rank() over(partition BY sample
                   ORDER BY ROW,col) AS ID,
       X,
       Y,
       H
FROM downsampled;

DELETE
FROM flowdata
WHERE H=0;

CREATE TABLE flowdata1 AS
SELECT S2.uniq_id,
       S1.sample,
       S1.id,
       S1.x,
       S1.y,
       S1.h
FROM
  (SELECT rank() over(
                      ORDER BY rank) AS Uniq_id,
                 rank
   FROM
     (SELECT rank
      FROM
        (SELECT *,
                rank() over(
                            ORDER BY X,Y)
         FROM flowdata) AS S
      GROUP BY rank) AS SS) AS S2,

  (SELECT *,
          rank() over(
                      ORDER BY X,Y)
   FROM flowdata) AS S1
WHERE S1.rank=S2.rank;

DROP TABLE flowdata;

ALTER TABLE flowdata1 RENAME TO flowdata;

CREATE TABLE
TEMPORARY AS
SELECT M.sample,
       M.neighbour,
       F.ID,
       F.uniq_id,
       F.X,
       F.Y
FROM macro_neighbour AS M
INNER JOIN flowdata AS F ON M.neighbour=F.sample;

/*
#<=100 denotes the distance along x and y axes. 200 is a preferred number fro the distance.
#rank<=150 indicates the number of neighbors you would like to store. It dictates the maximum size of the co-variance matrix.
*/
CREATE TABLE proximity AS
SELECT SAMPLE,
       SPATIALID,
       MACRO_NEIGHBOUR,
       SPATIAL_NEIGHBOUR,
       UNIQ_ID,
       UNIQ_ID_NEIGHBOUR
FROM
  (SELECT *,
          rank() over(partition BY sample,spatialid
                      ORDER BY random)
   FROM
     (SELECT F.sample,
             F.id AS spatialid,
             T.neighbour AS macro_neighbour,
             T.id AS spatial_neighbour,
             F.uniq_id,
             T.uniq_id AS uniq_id_neighbour,
             random() AS random
      FROM flowdata AS F,
      TEMPORARY AS T
      WHERE T.sample=F.sample
        AND abs(F.X-T.X)<=100
        AND abs(F.Y-T.Y)<=100) AS S) AS SS
WHERE rank<=150;

CREATE TABLE proximity_details AS
SELECT *
FROM
  (SELECT sample,
          spatialid,
          count(sample) over(partition BY sample,spatialid)
   FROM proximity) AS S
GROUP BY sample,
         spatialid,
         COUNT;

DROP TABLE
TEMPORARY;

CREATE TALBE resamples(resample int,P1 double PRECISION,P2 double PRECISION,P3 double PRECISION,P4 double PRECISION, w double PRECISION);
LOAD DATA LOCAL INFILE '/keith1/data/users/kmarcus2/sql-input/macro_resamples.tmp'
  INTO TABLE resamples
  FIELDS TERMINATED BY ' ';

CREATE TABLE maxmin_resamples AS
SELECT max(p1) AS maxp1,
       min(p1) AS minp1,
       max(p2) AS maxp2,
       min(p2) AS minp2,
       max(p3) AS maxp3,
       min(p3) AS minp3,
       max(p4) AS maxp4,
       min(p4) AS minp4
FROM resamples;

CREATE TABLE maxmin_samples AS
SELECT max(logvol) AS maxp1,
       min(logvol) AS minp1,
       max(direction) AS maxp2,
       min(direction) AS minp2,
       max(basal) AS maxp3,
       min(basal) AS minp3,
       max(internal) AS maxp4,
       min(internal) AS minp4
FROM meta;

CREATE TABLE maxmin AS
SELECT *
FROM maxmin_resamples
UNION
SELECT *
FROM maxmin_samples;

DROP TABLE maxmin_resamples,
           maxmin_samples;

CREATE TABLE minmax AS
SELECT min(minp1) AS minp1,
       max(maxp1) AS maxp1,
       min(minp2) AS minp2,
       max(maxp2) AS maxp2,
       min(minp3) AS minp3,
       max(maxp3) AS maxp3,
       min(minp4) AS minp4,
       max(maxp4) AS maxp4
FROM maxmin;

CREATE TABLE samples_scaled AS
SELECT sample,
       (S1.logvol-S2.minP1)/(maxP1 - minP1) AS P1,
       (S1.direction-S2.minP2)/(maxP2 - minP2) AS P2,
       (S1.basal-S2.minP3)/(maxP3 - minP3) AS P3,
       (S1.internal-S2.minP4)/(maxP4 - minP4) AS P4
FROM
  (SELECT sample,
          logvol,
          direction,
          basal,
          internal
   FROM meta) AS S1,

  (SELECT *
   FROM minmax) AS S2;

DROP TABLE maxmin;

CREATE TABLE resamples_scaled AS
SELECT resample,
       (S1.p1-S2.minP1)/(maxP1 - minP1) AS P1,
       (S1.p2-S2.minP2)/(maxP2 - minP2) AS P2,
       (S1.p3-S2.minP3)/(maxP3 - minP3) AS P3,
       (S1.p4-S2.minP4)/(maxP4 - minP4) AS P4
FROM
  (SELECT resample,
          P1,
          P2,
          P3,
          P4
   FROM resamples) AS S1,

  (SELECT *
   FROM minmax) AS S2;

CREATE TABLE Needed_resamples AS
SELECT sample,
       resample,
       sqrt((S1.P1-S2.P1)*(S1.P1-S2.P1) + (S1.P2-S2.P2)*(S1.P2-S2.P2) + (S1.P3-S2.P3)*(S1.P3-S2.P3) + (S1.P4-S2.P4)*(S1.P4-S2.P4)) AS distance
FROM
  (SELECT *
   FROM samples_scaled) AS S1,

  (SELECT *
   FROM resamples_scaled) AS S2
WHERE distance<.087;

CREATE TABLE new_scaled_resamples AS
SELECT SS.rank AS resample,
       S1.p1,
       S1.p2,
       S1.p3,
       S1.p4
FROM
  (SELECT resample AS old_resample,
          rank() over(
                      ORDER BY resample)
   FROM
     (SELECT resample
      FROM needed_resamples
      GROUP BY resample) AS S) AS SS,

  (SELECT *
   FROM resamples_scaled) AS S1
WHERE SS.old_resample=S1.resample;

CREATE TABLE new_resamples AS
SELECT SS.rank AS resample,
       S1.p1,
       S1.p2,
       S1.p3,
       S1.p4,
       S1.w
FROM
  (SELECT resample AS old_resample,
          rank() over(
                      ORDER BY resample)
   FROM
     (SELECT resample
      FROM needed_resamples
      GROUP BY resample) AS S) AS SS,

  (SELECT *
   FROM resamples) AS S1
WHERE SS.rank=S1.resample;

DROP TABLE resamples,
           resamples_scaled;


ALTER TABLE new_resamples RENAME TO resamples;

ALTER TABLE new_scaled_resamples RENAME TO resamples_scaled;

/* 0.12 dennotes the distance of search. You can change it to suit your need and your computational resource. */
CREATE TABLE resample_neighbours AS
SELECT *,
       (1/distance)/sum(1/distance) over(partition BY resample) AS weight
FROM
  (SELECT resample,
          sample,
          sqrt((S1.P1-S2.P1)*(S1.P1-S2.P1) + (S1.P2-S2.P2)*(S1.P2-S2.P2) + (S1.P3-S2.P3)*(S1.P3-S2.P3) + (S1.P4-S2.P4)*(S1.P4-S2.P4)) AS distance
   FROM
     (SELECT *
      FROM samples_scaled) AS S2,

     (SELECT *
      FROM resamples_scaled) AS S1
   WHERE distance<.12) AS S3;

DROP TABLE needed_resamples;

CREATE TABLE phm(id int,X double PRECISION, Y double PRECISION);
LOAD DATA LOCAL INFILE '/keith1/data/users/kmarcus2/sql-input/montserrat_take2_vol_dir_bed_int.phm'
  INTO TABLE phm
  FIELDS TERMINATED BY ' ';

CREATE TABLE uniq_coord AS
SELECT uniq_id,
       X,
       Y
FROM flowdata
GROUP BY uniq_id,
         X,
         Y;

CREATE TABLE
TEMPORARY AS
SELECT max(X) AS maxx,
       min(x) AS minx,
       max(y) AS maxy,
       min(y) AS miny
FROM uniq_coord;

CREATE TABLE phm_needed AS
SELECT P.id,
       P.x,
       P.y
FROM phm AS P,
TEMPORARY AS T
WHERE P.X <= T.maxx
  AND P.X >= T.minx
  AND P.Y <= T.maxy
  AND P.Y >= T.miny;

DROP TABLE
TEMPORARY;


CREATE TABLE phm_neighbours AS
SELECT *
FROM
  (SELECT *,
          rank() over(partition BY sample,phm_id
                      ORDER BY distance)
   FROM
     (SELECT S2.sample,
             S2.uniq_id,
             S2.phm_id,
             S4.spatialid,
             S2.distance
      FROM
        (SELECT F.sample,
                F.id AS spatialid,
                F.uniq_id,
                S1.phm_id,
                S1.distance
         FROM
           (SELECT P.id AS phm_id,
                   U.uniq_id,
                   sqrt((P.X-U.X)*(P.X-U.X)+(P.Y-U.Y)*(P.Y-U.Y)) AS distance
            FROM uniq_coord AS U,
                 phm_needed AS P
            WHERE distance<=100) AS S1
         INNER JOIN flowdata AS F ON F.uniq_id=S1.uniq_id) AS S2,

        (SELECT S3.sample,
                S3.spatialid,
                F.uniq_id,
                S3.indicator
         FROM
           (SELECT *,
                   COUNT/10 AS indicator
            FROM proximity_details) AS S3,
              flowdata AS F
         WHERE S3.sample=F.sample
           AND S3.spatialid=F.id
           AND S3.indicator>0) AS S4
      WHERE S2.sample=S4.sample
        AND S2.uniq_id=S4.uniq_id) AS S) AS SS
WHERE SS.rank<=3;

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

CREATE TABLE sample_count AS
SELECT sample,
       sum(COUNT) over(
                       ORDER BY sample ASC ROWS unbounded preceding exclude CURRENT ROW)
FROM
  (SELECT *
   FROM
     (SELECT sample,
             COUNT(sample) over(partition BY sample)
      FROM proximity) AS S
   GROUP BY sample,
            COUNT) AS SS;

UPDATE sample_count
SET SUM=0
WHERE SUM IS NULL;

CREATE TABLE spatial_count AS
SELECT *,
       SUM(COUNT) over(partition BY sample
                       ORDER BY spatialid ASC ROWS unbounded preceding exclude CURRENT ROW)
FROM
  (SELECT SS.sample,
          SS.spatialid,
          F.H,
          SS.COUNT
   FROM
     (SELECT *
      FROM
        (SELECT sample,
                spatialid,
                COUNT(spatialid) over(partition BY sample,spatialid)
         FROM proximity) AS S
      GROUP BY S.spatialid,
               S.sample,
               S.COUNT) AS SS
   INNER JOIN flowdata AS F ON SS.sample=F.sample
   AND SS.spatialid=F.id
   ORDER BY sample,
            spatialid) AS SSS;

UPDATE spatial_count
SET SUM=0
WHERE SUM IS NULL;

CREATE TABLE weighted_phm(sample int, resample int, phm_id int, mean double PRECISION, variance double PRECISION) distribute ON (phm_id,
                                                                                                                                 resample,
                                                                                                                                 sample);
