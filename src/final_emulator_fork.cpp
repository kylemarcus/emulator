#include "mpi.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <vector>
#include <list>
#include <unistd.h>
#include <sys/wait.h>
#include "build_emulator.h"
#include "matrixop.h"
#include <ctime>
#include <string.h>

#define RESAMPLES_TXT "/panasas/scratch/kmarcus2/emulator/extract/resamples.txt"
#define PHM_TXT "/panasas/scratch/kmarcus2/emulator/extract/phm.txt"
#define UNCERTAIN_INPUT_LIST_TXT "/panasas/scratch/kmarcus2/emulator/extract/uncertain_input_list.txt"
#define SAMPLE_COUNT_TXT "/panasas/scratch/kmarcus2/emulator/extract/sample_count.txt"
#define UNIQ_COORDS_TXT "/panasas/scratch/kmarcus2/emulator/extract/uniq_coords.txt"

#define SPATIAL_COUNT_SH "/user/kmarcus2/emulator/src/scripts/spatial_count.sh"
#define CREATE_TABLE_SH "/user/kmarcus2/emulator/src/scripts/create_table.sh"
#define NZLOAD_SH "/user/kmarcus2/emulator/src/scripts/nzload.sh"
#define PHM_RESAMPLE_SH "/user/kmarcus2/emulator/src/scripts/phm_resample.sh"
#define PROXIMITY_SH "/user/kmarcus2/emulator/src/scripts/proximity.sh"
#define STATS_SH "/user/kmarcus2/emulator/src/scripts/stats.sh"

#define RADIUS 100						//Radius of neighbour search for spatial points
#define Ndim 4							//Number of dimensions - parametric
#define M 32							//Number of processors in a group
#define S 2048							//Number of Titan simulations
#define MACRO_SCALED_RADIUS 0.2			//Radius of search for macro-parameters
#define T 120							//Number of neighbours needed to be considered for setting up co-variance matrix. should be < 500
#define MAX_COUNT 800					//Number of emulator constructions after which data is piped to Netezza
#define TOO_FEW 10						//Number of neighbours for a spatial point. for 10 points or less 
#define MERGE_COUNT 200
#define MERGE_SAFE 25
#define NEWTON_ITERATIONS 10

using namespace std;

struct Data6
{
int Sample,Uniq_record,Return_count;
double *mean_and_variance;
};

void BuildX(double * X, struct Data4 * proximity, vector<struct Data3> & uniq_coord,struct Data5 & I,int i,struct Data2 * q,struct Data1 *p,int sample);
void Send_To_Netezza(list<struct Data6> & Pipeout, int rank,int group, int fid);

int main(int argc, char * argv[])
{

int rank,Nprocs;
MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD,&rank);
MPI_Comm_size(MPI_COMM_WORLD,&Nprocs);
MPI_Status root_stat;
MPI_Request root_req;
int group=0,root,start_rank,stop_rank,Np;
int i,j,k,m,n,sample=0,trash,skip,indicator;
double garb,lmin[Ndim];
char str[200],filename[50],merge_file[50];
double *min,*max;
FILE * fp1, * fpipe;
double **resamples,**phm;
int resample_neighbours_count, * resample_neighbours, * phm_neighbours_count, **phm_neighbours,*phm_neighbours_index; //phm_neighbours_index is for rank != root while phm_neighbours is for rank == root
int Ngroups,Nleft,Nsample,Ntasks,start_sample,stop_sample,*rank_list;
int isample,iproc,Netezza_count=0;
int count, **pass_count;

struct Data5 I;
struct Data1 *p;
struct Data2 *q;
struct Data4 *r;
struct Data6 **Netezza;
list<struct Data6> Pipeout;
list<struct Data6> :: iterator pipe_it;

Ngroups = int(Nprocs/M);
Nleft = Nprocs-Ngroups*M;		//----------------------------------------------------> MODIFIED ON JULY 16 ON 1:31 pm
cout << " NPROCS = " << Nprocs << "\n";
//ONLY FOR RANK = 0
int master_group[Ngroups],master_rank[Ngroups],master_flag[Ngroups],master_indicator[Ngroups],master_sum=S,prev_sample[Ngroups];
int Netz[Ngroups];
MPI_Request master_req[Ngroups];
MPI_Status master_stat[Ngroups];

for(i=1;i<=Ngroups;i++)
Netz[i] = -1;

for(i=1;i<=Ngroups;i++)
{
	m = i <= Nleft ? i-1 : Nleft;
	n = i <= Nleft ? 0 : (i-Nleft-1);
	k = i <= Nleft ? M+1 : M;
	Ntasks = k;

	if ( rank >= m*(M+1)+n*M && rank <= m*(M+1)+n*M + k )			//////////////////// MODIFIED RIGHT NOW JULY 16 1:31 pm ///////////////////////////
	{
		root = m*(M+1)+n*M;               // ----------------------------------------------------> MODIFIED ON JULY ON 1:31 pm
		start_rank = m*(M+1)+n*M + 1;	  // ----------------------------------------------------> MODIFIED ON JULY ON 1:31 pm
		stop_rank = m*(M+1)+n*M + k - 1;  // ----------------------------------------------------> MODIFIED ON JULY ON 1:31 pm
		group = i;

	}

	if ( rank == 0 )
	{
		master_rank[i-1] = m*(M+1)+n*M > 0 ? m*(M+1)+n*M : m*(M+1)+n*M + 1;
		master_group[i] = i;
	}
}


if ( group == 1 )
{
root = 1;
start_rank = 2;
Ntasks = Ntasks -1;
}

//master_sum = 1;
///////////////////////////////////////////////////    ADDED ON JULY 16 -- 1:36 pm   ///////////////////////////////////////////////////////
if ( rank == 0 )
{
bool merge[S];
int Mi,Mv=1,Msum=0;

	for(j=1;j<S;j++)
	merge[j] = 0;

	sample = 0;
	for(j=1;j<=Ngroups;j++)
        {
		i = j-1;
                MPI_Send(&sample,1,MPI_INT,master_rank[i],8,MPI_COMM_WORLD);
                sample++;

		MPI_Send(&sample,1,MPI_INT,master_rank[i],7,MPI_COMM_WORLD);            ///////////////////    ADDED on JULY 17 --- 1:25 pm //////////////
                sample++;

                master_flag[i] = 0;
		master_indicator[i] = -1;
                MPI_Irecv(&master_indicator[i],1,MPI_INT,master_rank[i],Nprocs+master_rank[i],MPI_COMM_WORLD,&master_req[i]);
        }

while ( master_sum > 0 )
{
	for(j=1;j<=Ngroups;j++)
	{
		i = j-1;
	        MPI_Test(&master_req[i],&master_flag[i],&master_stat[i]);

	        if ( master_flag[i] == 1 && master_indicator[i] > -1 )
	        {

                        master_indicator[i] = -1;
	                master_sum--;
	                MPI_Send(&sample,1,MPI_INT,master_rank[i],7,MPI_COMM_WORLD);
	                sample++;

	                master_flag[i] = 0;
	                MPI_Irecv(&master_indicator[i],1,MPI_INT,master_rank[i],Nprocs+master_rank[i],MPI_COMM_WORLD,&master_req[i]);

	        }


/*
		if ( master_flag[i] == 1 && master_indicator[i] == -2 )
		{
			sprintf(stat_file,"./stats.sh %d &",j);
			system(stat_file);
	                MPI_Irecv(&master_indicator[i],1,MPI_INT,master_rank[i],Nprocs+master_rank[i],MPI_COMM_WORLD,&master_req[i]);
			master_flag[i] = 0;
			master_indicator[i] = -1;
		}
*/

	}
}


//for(i=1;i<=Ngroups;i++)
//MPI_Send(&sample,1,MPI_INT,master_rank[i-1],7,MPI_COMM_WORLD);

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


Np = Ntasks;
Ntasks = Ntasks - 1;

k = S - int(S/Ngroups)*Ngroups;
Nsample = group <= k ? int(S/Ngroups)+1 : int(S/Ngroups);

j = 0;
if ( rank == root && rank != 0 )
{
	char create_table[25];
	sprintf(create_table,CREATE_TABLE_SH" %d",group);
	system(create_table);

	Netezza = new struct Data6*[Ntasks];

	pass_count = new int * [Ntasks];
	for(i=0;i<Ntasks;i++)
	pass_count[i] = new int [3];

	MPI_Recv(&sample,1,MPI_INT,0,8,MPI_COMM_WORLD,&root_stat);          /////////////////////////   ADDED ON JULY 16 ---- 1:15 pm ///////////////////////
	cout << " sample = " << sample;
        start_sample = sample;
        MPI_Irecv(&isample,1,MPI_INT,0,7,MPI_COMM_WORLD,&root_req);

	rank_list = new int[Ntasks];
	for(i=start_rank;i<=stop_rank;i++)
	rank_list[j++] = i;
}
else
rank_list = NULL;


if ( rank != root && rank != 0 )
{
pass_count = new int *;
pass_count[0] = new int[3];

///////////////////////////////// Read resamples.txt and phm.txt and store the data //////////////////////////////////
///////////////////////////////// all computing processors store their copy  /////////////////////////////////////////
fp1 = fopen(RESAMPLES_TXT,"r");
fscanf(fp1,"%d",&k);

resamples = new double * [k];
for(i=0;i<k;i++)
resamples[i] = new double[Ndim+2];

for(i=0;i<k;i++)
{
	fscanf(fp1,"%lf,",&garb);
	for(j=0;j<Ndim;j++)
	fscanf(fp1,"%lf,",&resamples[i][j]);
	fscanf(fp1,"%lf",&garb);
}
fclose(fp1);

fp1 = fopen(PHM_TXT,"r");
fscanf(fp1,"%d",&k);

phm = new double * [k];
for(i=0;i<k;i++)
phm[i] = new double[2];

for(i=0;i<k;i++)
fscanf(fp1,"%lf,%lf",&phm[i][0],&phm[i][1]);

fclose(fp1);
}			//end of if( rank != root )


if ( rank == root && rank != 0 )
{

MPI_Request req1[Ntasks],req2[Ntasks],req3[Ntasks],req4[Ntasks],req5[Ntasks];
MPI_Status stat1[Ntasks],stat2[Ntasks],stat3[Ntasks],stat4[Ntasks],stat5[Ntasks];

int fork_count=0,fid;
pid_t fid_fork;
int stat_count=0;
double **X;
int icount=0,fifo;
int sample_record[Ntasks],uniq_record[Ntasks],return_count[Ntasks],index;
int flag[Ntasks];
struct Data4 * proximity;
X = new double * [Ntasks];

char stat_file[25],Loading_file[25];

//////////////////////////////////////////    ADDED ON JULY 22 2012 ----------- 5:37pm   ///////////////////////////////////
char pipename[25];

sprintf(pipename,"Netpipe%d",rank);
//sprintf(Loading_file,"./nzload.sh %d %d",rank,group);

fifo = mkfifo(pipename,0666);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////// Read uncertain_input_list.txt and store it  ///////////////////////////////////////////////////
fp1 = fopen(UNCERTAIN_INPUT_LIST_TXT,"r");
fscanf(fp1,"additional file format lines=%d\n",&skip);
for(i=0;i<skip;i++)
fgets(str,200,fp1);
//fscanf(fp1,"%d",&Ndim);
fscanf(fp1,"%d",&trash);
fscanf(fp1,"%d",&trash);

I.data = new double *[2048];
for(i=0;i< 2048;i++)
I.data[i] = new double [Ndim];

min = new double [2*Ndim];
max = new double [2*Ndim];

for(i=0;i<2048;i++)
{
	for(j=0;j<Ndim;j++)
	{
		fscanf(fp1,"%lf ",&(I.data[i][j]));

		if ( i == 0 )
		{
			min[j] = I.data[i][j];
			max[j] = I.data[i][j];
		}

		else
		{
			if ( I.data[i][j] < min[j] )
			min[j] = I.data[i][j];
			if ( I.data[i][j] > max[j] )
			max[j] = I.data[i][j];
		}
	}
}
fclose(fp1);

for(i=0;i<Ndim;i++)
lmin[i] = 0.25*MACRO_SCALED_RADIUS*(max[i]-min[i]);

for(i=0;i<Ntasks;i++)
MPI_Send(lmin,4,MPI_DOUBLE,rank_list[i],5,MPI_COMM_WORLD);

delete [] max;
delete [] min;

////////////////////////   Read sample_count.txt and spatial_count.txt and store it in MetaData structure ////////////////////////////
fp1 = fopen(SAMPLE_COUNT_TXT,"r");
fpipe = popen(SPATIAL_COUNT_SH,"r");

struct Data1 * MetaData = new struct Data1[S];

for(i=0;i<S;i++)
{
	fscanf(fp1,"%d,%d",&trash,&MetaData[i].spatial_count);

	MetaData[i].spatial = new struct Data2[MetaData[i].spatial_count];
	p = &MetaData[i];

	for(j=0;j< p->spatial_count;j++)
	{
		q = &(p->spatial[j]);
		fscanf(fpipe,"%lf,%d", &(q->H), &(q->neighbour_count) );

	}
}

fclose(fp1);
pclose(fpipe);

///////////////////////////////////   Read and store uniq_coords   ///////////////////////////////////////////////////////
struct Data3 temp;
vector<struct Data3> uniq_coord;
fp1 = fopen(UNIQ_COORDS_TXT,"r");

while( fscanf(fp1,"%d,%lf,%lf",&trash,&(temp.x),&(temp.y) ) != EOF )
uniq_coord.push_back(temp);

fclose(fp1);

bool track[Ntasks];

///////////////////// FROM HERE EVERYTHing SHOULD BE LOOPED////////////////////////////////

while ( sample < S )
{
p = &MetaData[sample];
q = p->spatial;

if ( sample != start_sample )
{
MPI_Waitall(Ntasks,req1,stat1);
MPI_Waitall(Ntasks,req2,stat2);
MPI_Waitall(Ntasks,req3,stat3);
MPI_Waitall(Ntasks,req4,stat4);
MPI_Waitall(Ntasks,req5,stat5);
}

//////////////////////// Read from netezza and create data structures for storing resample_neighbours and phm_neighbours  ////////////////////
sprintf(filename,PHM_RESAMPLE_SH" %d",sample+1);
fpipe = popen(filename,"r");
fscanf(fpipe,"%d",&resample_neighbours_count);
resample_neighbours = new int [resample_neighbours_count];

for(i=0;i<resample_neighbours_count;i++)
fscanf(fpipe,"%d",&resample_neighbours[i]);

phm_neighbours_count = new int [p->spatial_count];
for(i=0; i< p->spatial_count; i++)
fscanf(fpipe,"%d",&phm_neighbours_count[i]);

phm_neighbours = new int * [p->spatial_count];

for(i=0;i< p->spatial_count;i++)
{
	phm_neighbours[i] = new int [phm_neighbours_count[i]];

	for(j=0;j< phm_neighbours_count[i]; j++)
	fscanf(fpipe,"%d",&phm_neighbours[i][j]);

}

pclose(fpipe);


//////////////////////// Read and store proximity information for the sample number in progress //////////////////////////////////////////////
sprintf(filename,PROXIMITY_SH" %d",sample+1);			//sample+1 because C uses index from 0 and Netezza uses index from 1. The variable
cout << filename << "\n";
fpipe = popen(filename,"r");					//sample in this code refers to the index in C and so starts from 0
proximity = new struct Data4[MetaData[sample].spatial_count];
p = &(MetaData[sample]);
q = p->spatial;
r = proximity;

for(i=0; i< p->spatial_count ; i++)
{

	r->spatial_neighbour = new int[q->neighbour_count];
	r->uniq_id_neighbour = new int[q->neighbour_count];
	r->macro_neighbour = new int[q->neighbour_count];

	for(j=0;j< q->neighbour_count;j++)
	fscanf(fpipe,"%d,%d,%d,%d",&(r->uniq_id),&(r->spatial_neighbour[j]),&(r->uniq_id_neighbour[j]),&(r->macro_neighbour[j]));

	q++;
	r++;

}

pclose(fpipe);

//Pointers are pointed afresh to start the loop over spatial points --- Loop over spatial points begins here
p = &(MetaData[sample]);
q = p->spatial;
r = proximity;
icount = 0;

if ( sample == start_sample )
{

cout << " rank = " << rank << " I AM HERE " << __LINE__ << "\n";
for(iproc=0;iproc<Ntasks;iproc++)
{

/////////////////// ADDED NOW JULY13 ---- 2:17 pm //////////////////////////////
	while ( (phm_neighbours_count[icount] == 0) || (q->neighbour_count < TOO_FEW) )
	{
		q++;
		icount++;

	}

///////////////////////////////////////////////////////////////////////////////////

	indicator = 2;
	MPI_Isend(&indicator,1,MPI_INT,rank_list[iproc],0,MPI_COMM_WORLD,&req1[iproc]);

	Netezza[iproc] = new struct Data6;
	Netezza[iproc][0].Sample  = sample;
	Netezza[iproc][0].Uniq_record = proximity[icount].uniq_id;
	Netezza[iproc][0].Return_count = 4*resample_neighbours_count*phm_neighbours_count[icount];

	sample_record[iproc] = sample;
	uniq_record[iproc] = proximity[icount].uniq_id;
	return_count[iproc] = resample_neighbours_count*phm_neighbours_count[icount];

	count = q->neighbour_count;
	X[iproc] = new double[(Ndim+3)*(count+1)];
	BuildX(&X[iproc][0],proximity,uniq_coord,I,icount,q,p,sample);

	pass_count[iproc][0] = resample_neighbours_count;
	pass_count[iproc][1] = phm_neighbours_count[icount];
	pass_count[iproc][2] = q->neighbour_count;

	MPI_Isend(&pass_count[iproc][0],3,MPI_INT,rank_list[iproc],3,MPI_COMM_WORLD,&req2[iproc]);

	MPI_Isend(resample_neighbours,resample_neighbours_count,MPI_INT,rank_list[iproc],1,MPI_COMM_WORLD,&req3[iproc]);

	MPI_Isend(&phm_neighbours[icount][0],phm_neighbours_count[icount],MPI_INT,rank_list[iproc],2,MPI_COMM_WORLD,&req4[iproc]);

	MPI_Isend(&X[iproc][0],(q->neighbour_count)*(Ndim+3),MPI_DOUBLE,rank_list[iproc],4,MPI_COMM_WORLD,&req5[iproc]);

	Netezza[iproc][0].mean_and_variance = new double[4*return_count[iproc]];
	MPI_Irecv(Netezza[iproc][0].mean_and_variance,4*return_count[iproc],MPI_DOUBLE,rank_list[iproc],Nprocs+rank_list[iproc],MPI_COMM_WORLD,&req2[iproc]);

//	mean_and_variance[iproc] = new double[4*return_count[iproc]];
//	MPI_Irecv(&mean_and_variance[iproc][0],4*return_count[iproc],MPI_DOUBLE,rank_list[iproc],Nprocs+rank_list[iproc],MPI_COMM_WORLD,&req2[iproc]);

	q++;
	icount++;


}

}


while ( icount < p->spatial_count )
{

  //if ( sample ==7 )
  //cout << p->spatial_count << "\n";

	for(iproc=0;iproc<Ntasks;iproc++)
	{
		MPI_Test(&req2[iproc],&flag[iproc],&stat2[iproc]);

/////////////////// ADDED NOW JULY13 ---- 2:17 pm //////////////////////////////
		while ( (phm_neighbours_count[icount] == 0) || (q->neighbour_count < 10) )
		{

			q++;
			icount++;

			if ( icount == p->spatial_count )
			break;

		}

		if ( icount >= p->spatial_count )
		break;
///////////////////////////////////////////////////////////////////////////////////


		if ( flag[iproc] == 1 )
		{
			track[iproc] = 1;

			if ( sample_record[iproc] != sample )
			indicator = 2;
			else
			indicator = 1;

			MPI_Isend(&indicator,1,MPI_INT,rank_list[iproc],0,MPI_COMM_WORLD,&req1[iproc]);
//			cout << "  rank = " << rank << "   I  AM HERE " << __LINE__ << "   icount = " << icount << "     sample = " << sample << "     Netezza_count = " << Netezza_count <<"\n";

			delete [] X[iproc];
			count = q->neighbour_count;

			X[iproc] = new double[(Ndim+3)*(count+1)];
			BuildX(&X[iproc][0],proximity,uniq_coord,I,icount,q,p,sample);

			pass_count[iproc][0] = resample_neighbours_count;
			pass_count[iproc][1] = phm_neighbours_count[icount];
			pass_count[iproc][2] = q->neighbour_count;

			MPI_Isend(&pass_count[iproc][0],3,MPI_INT,rank_list[iproc],3,MPI_COMM_WORLD,&req2[iproc]);

			if ( indicator == 2 )
			MPI_Isend(resample_neighbours,resample_neighbours_count,MPI_INT,rank_list[iproc],1,MPI_COMM_WORLD,&req3[iproc]);

			MPI_Isend(&phm_neighbours[icount][0],phm_neighbours_count[icount],MPI_INT,rank_list[iproc],2,MPI_COMM_WORLD,&req4[iproc]);

			MPI_Isend(&X[iproc][0],(q->neighbour_count)*(Ndim+3),MPI_DOUBLE,rank_list[iproc],4,MPI_COMM_WORLD,&req5[iproc]);

			Pipeout.push_back(Netezza[iproc][0]);
			Netezza_count++;

			if ( Netezza_count == MAX_COUNT )
			{

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// SEPT 5, 2012 COMMENTING BELOW LINES AND ADDING MODIFIED PARTS OF BELOW LINES   ///////////////
///////////////////////// BECAUSE 8 CORE NODES DISPLA COME HYBRID ERROR DUE TO SPAWNING OF CHILD PROCES    //////////////
/*
				cout << " fork_fid = " << fid_fork << "   fork_count =" << fork_count << "\n";
				if ( fork_count > 0 && fid_fork > 0 )
				{

					wait(NULL);
//					MPI_Send(&k,1,MPI_INT,0,Nprocs+rank,MPI_COMM_WORLD);
					cout << " FINISHED WAITING .....\n";
				}


				if ( stat_count == 10 )
				{
					close(fid);
					sprintf(stat_file,"./stats.sh %d",group);
					system(stat_file);
					stat_count = 0;
				}
*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

				if ( stat_count == 10 )
				{
					close(fid);
					sprintf(stat_file,STATS_SH" %d",group);
					system(stat_file);
					stat_count = 0;
				}

				stat_count++;

				if ( stat_count == 1 )
				{
					sprintf(Loading_file,NZLOAD_SH" %d %d",rank,group);
					sprintf(pipename,"Netpipe%d",rank);
					system(Loading_file);
					fid = open(pipename,O_WRONLY);

					if ( fid < 0 )
					cout << " COULD NOT OPEN THE PIPE FOR WRITING \n";

				}

////////////////////////////////////  SEPT 5,2012 NEW LINES ADDED BELOW-- DO NOT DELETE //////////////////////////////////////////////

				Send_To_Netezza(Pipeout,rank,group,fid);
                                Netezza_count = 0;
                                Pipeout.erase(Pipeout.begin(),Pipeout.end());

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			}   //end of if( Netezza_count == MAX_COUNT )


			Netezza[iproc] = new struct Data6;
			Netezza[iproc][0].Sample  = sample;
			Netezza[iproc][0].Uniq_record = proximity[icount].uniq_id;
			Netezza[iproc][0].Return_count = 4*resample_neighbours_count*phm_neighbours_count[icount];

			sample_record[iproc] = sample;
			uniq_record[iproc] = proximity[icount].uniq_id;
			return_count[iproc] = resample_neighbours_count*phm_neighbours_count[icount];

			Netezza[iproc][0].mean_and_variance = new double[4*return_count[iproc]];
			MPI_Irecv(Netezza[iproc][0].mean_and_variance,4*return_count[iproc],MPI_DOUBLE,rank_list[iproc],Nprocs+rank_list[iproc],MPI_COMM_WORLD,&req2[iproc]);

//			delete [] mean_and_variance[iproc];
//			mean_and_variance[iproc] = new double[4*return_count[iproc]];
//			MPI_Irecv(&mean_and_variance[iproc][0],4*return_count[iproc],MPI_DOUBLE,rank_list[iproc],Nprocs+rank_list[iproc],MPI_COMM_WORLD,&req2[iproc]);

			track[iproc] = 1;

			q++;
			icount++;
			flag[iproc] = 0;

//SEPT 17,2012 The two lines below are repeated to skip emulator construction on every alternate point. I'm doing this because I'm increasing
//the number of downsampled points. The idea is to have sufficients points in neighbourhood but since more points lead to more
// emulator constructions which is uncecessary every alternate point will be skipped. Currently 100 metres search radius results in no spatial
// point other than the point itself. This leads to rank-deficient polynomial matrix G which is obvious. It is not possible to do linear regression
// if there is only one point in any of the dimensions.

//Nov12, 2012 Skipping the points twice because downsampling has selected close to 6000 points per sample.

//			if ( p->spatial_count > 2500 )
//			{
				q++;
				icount++;
				q++;
                                icount++;

//			}

////////////////////////////////////////////////////////////////////////////////////////////

/*			if ( icount == 80 )
			icount = p->spatial_count;
*/
///////////////////////////////////////////////////////////////////////////////////////////////

		}

		if ( icount == p->spatial_count )
		break;

	}	//end of for loop

}		//end of while loop



for(i=0;i < p->spatial_count;i++)
{
	delete [] phm_neighbours[i];
	delete [] proximity[i].spatial_neighbour;
	delete [] proximity[i].uniq_id_neighbour;
	delete [] proximity[i].macro_neighbour;
}
delete [] resample_neighbours;
delete [] phm_neighbours;
delete [] phm_neighbours_count;
delete [] proximity;

/////////////////////////   ADDED ON JULY 16 ---- 1:15 pm ///////////////////////
MPI_Send(&sample,1,MPI_INT,0,Nprocs+rank,MPI_COMM_WORLD);
sample = isample;
if ( sample >= S )
continue;

MPI_Wait(&root_req,&root_stat);
MPI_Irecv(&isample,1,MPI_INT,0,7,MPI_COMM_WORLD,&root_req);

/////////////////////////////////////////////////////////////////////////////////

}			//end of sample loop


int sum = 0;

for(iproc=0;iproc < Ntasks;iproc++)
sum = track[iproc] + sum;
//sum = flag[iproc] + sum;
indicator = 0;

while ( sum != 0 )
{

	sum = 0;
	for(iproc=0;iproc < Ntasks ;iproc++)
	{
		if ( track[iproc] == 1 )
		MPI_Test(&req2[iproc],&flag[iproc],&stat2[iproc]);

		if ( flag[iproc] == 1 && track[iproc] == 1 )
		{

			MPI_Send(&indicator,1,MPI_INT,rank_list[iproc],0,MPI_COMM_WORLD);

			track[iproc] = 0;

			Pipeout.push_back(Netezza[iproc][0]);
			Netezza_count++;

			flag[iproc] = 0;
//			delete [] mean_and_variance[iproc];
			delete [] X[iproc];

		}

		sum = sum + track[iproc];
	}

}

/*
fid_fork = fork();

if ( fid_fork == 0 )
{
	Send_To_Netezza(Pipeout,rank,group,fid);
	Pipeout.erase(Pipeout.begin(),Pipeout.end());			//AUGUST 25
	exit(1);
}

if ( fid_fork > 0 )
{
	wait(NULL);

	close(fid);
	sprintf(stat_file,"./stats.sh %d",group);
	system(stat_file);

	pipe_it  = Pipeout.begin();
	while( pipe_it != Pipeout.end() )
	{
		delete [] pipe_it->mean_and_variance;
		pipe_it++;
	}

	Pipeout.erase(Pipeout.begin(),Pipeout.end());
}
*/

//if ( fid <= 0 )
//fid = open(pipename,O_WRONLY);


//fid = open(pipename,O_WRONLY);

close(fid);

sprintf(Loading_file,NZLOAD_SH" %d %d",rank,group);
sprintf(pipename,"Netpipe%d",rank);
system(Loading_file);

fid = open(pipename,O_WRONLY);
Send_To_Netezza(Pipeout,rank,group,fid);
Pipeout.erase(Pipeout.begin(),Pipeout.end());
close(fid);

sprintf(stat_file,STATS_SH" %d",group);
system(stat_file);


for(iproc=0;iproc<Ntasks;iproc++)
delete Netezza[iproc];
delete [] Netezza;

//MPI_Irecv(&isample,1,MPI_INT,0,7,MPI_COMM_WORLD,&root_req);

}		//end of if ( root == rank )



if ( rank != root && rank != 0 )
{
MPI_Status stat[6],stat_wait,stat_ind;
MPI_Request req_wait,req_ind;
vector<double> mean_and_variance_vector;
double * mean_var;
struct Data5 Xsend,Ysend;

double *Y,*X,sigma,temp_resample[Ndim];
double log_min_corr_len[6],log_max_corr_len[6],log_corr_len_guess[6];
matrixop R,beta,rhs,G,Rinv;
iproc = 0;

MPI_Recv(lmin,Ndim,MPI_DOUBLE,root,5,MPI_COMM_WORLD,&stat[5]);
resample_neighbours = NULL;

log_min_corr_len[0] = log(RADIUS/4);
log_min_corr_len[1] = log(RADIUS/4);
log_max_corr_len[0] = log(2*RADIUS);
log_max_corr_len[1] = log(2*RADIUS);

for(i=0;i<Ndim;i++)
{
	log_min_corr_len[i+2] = log(lmin[i]);
	log_max_corr_len[i+2] = log(8*lmin[i]);
}

	MPI_Recv(&indicator,1,MPI_INT,root,0,MPI_COMM_WORLD,&stat[0]);

	while ( indicator != 0 )
	{
		MPI_Recv(&pass_count[0][0],3,MPI_INT,root,3,MPI_COMM_WORLD,&stat[3]);

		if ( indicator == 2 )
		{
			if ( resample_neighbours != NULL )
			delete [] resample_neighbours;

			resample_neighbours = new int [pass_count[0][0]];
			MPI_Recv(resample_neighbours,pass_count[0][0],MPI_INT,root,1,MPI_COMM_WORLD,&stat[1]);

		}

		phm_neighbours_index = new int[pass_count[0][1]];
		X = new double [pass_count[0][2]*(Ndim+3)];

		MPI_Recv(phm_neighbours_index,pass_count[0][1],MPI_INT,root,2,MPI_COMM_WORLD,&stat[2]);

		MPI_Recv(X,pass_count[0][2]*(Ndim+3),MPI_DOUBLE,root,4,MPI_COMM_WORLD,&stat[4]);

		count = pass_count[0][2] <= T ? pass_count[0][2]:T;
		Ysend.data = new double * [count];
		for(i=0;i<count;i++)
		Ysend.data[i] = new double[Ndim+2];
		Y = new double[count];
		k = 0;

		for(i=0;i<count;i++)
		{
			for(j=0;j<Ndim+2;j++)
			Ysend.data[i][j] = X[k++];
			Y[i] = X[k++];
		}

		log_corr_len_guess[0] = log(RADIUS);
		log_corr_len_guess[1] = log(RADIUS);

		for(i=0;i<Ndim;i++)
		log_corr_len_guess[i+2] = log(4*lmin[i]);

		for (i=0;i<NEWTON_ITERATIONS;i++)
		sigma = newton_method(R,Rinv,G,beta,rhs,i,Ndim+2,count,Ysend,Y,log_corr_len_guess,log_max_corr_len,log_min_corr_len,rank);
//		sigma = newton_method(R,Rinv,G,beta,rhs,i,Ndim+2,pass_count[0][2],Ysend,Y,log_corr_len_guess,log_max_corr_len,log_min_corr_len,rank);


		Xsend.data = new double * [pass_count[0][1]];
		for(i=0;i<pass_count[0][1];i++)
		Xsend.data[i] = new double [Ndim+2];

		double ymax=0;
		for(i=0;i<count;i++)
		{
			if ( ymax < Y[i] )
			ymax = Y[i];
		}

		for(i=0;i<pass_count[0][0];i++)			//loop over resamples -- i.e.calculate mean and variance using phm_neighbours for every resample neighbour
		{
			for(j=0;j<Ndim;j++)
			temp_resample[j] = resamples[resample_neighbours[i]][j];



			for(j=0;j<pass_count[0][1];j++)
			{
//				cout << phm_neighbours_index[j] << "   ";

				Xsend.data[j][0] = phm[phm_neighbours_index[j]][0];
				Xsend.data[j][1] = phm[phm_neighbours_index[j]][1];

				for(k=2;k<Ndim+2;k++)
				Xsend.data[j][k] = temp_resample[k-2];
			}

			Mean_and_Variance(mean_and_variance_vector, G, rhs, R, Rinv, beta, Xsend, Ysend, pass_count[0][1], count, Ndim+2, sigma, log_corr_len_guess, ymax);

		}

		if ( iproc > 0 )
		{
			MPI_Wait(&req_wait,&stat_wait);
			delete [] mean_var;
		}

		mean_var = new double [2*mean_and_variance_vector.size()];

		k = 0;
		m = 0;

		for(i=0;i<pass_count[0][0];i++)
		{
			for(j=0;j<pass_count[0][1];j++)
			{
				if ( mean_and_variance_vector[m] != -1 )
				{
					mean_var[k++] = resample_neighbours[i]+1;
					mean_var[k++] = phm_neighbours_index[j]+1;
					mean_var[k++] = mean_and_variance_vector[m++];
					mean_var[k++] = mean_and_variance_vector[m++];
				}

				else
				{
					mean_var[k++] = -1;
					mean_var[k++] = -1;
					mean_var[k++] = -1;
					mean_var[k++] = -1;
				}

			}
		}

//		MPI_Send(mean_var,k,MPI_DOUBLE,root,Nprocs+rank,MPI_COMM_WORLD);
		MPI_Isend(mean_var,k,MPI_DOUBLE,root,Nprocs+rank,MPI_COMM_WORLD,&req_wait);

		mean_and_variance_vector.erase(mean_and_variance_vector.begin(),mean_and_variance_vector.end());

//		MPI_Recv(&indicator,1,MPI_INT,root,0,MPI_COMM_WORLD,&stat[0]);


		MPI_Irecv(&indicator,1,MPI_INT,root,0,MPI_COMM_WORLD,&req_ind);
//		if ( iproc > 0 )
		MPI_Wait(&req_ind,&stat_ind);

		iproc++;

		for(i=0;i<pass_count[0][1];i++)
		delete [] Xsend.data[i];
		delete [] Xsend.data;
		for(i=0;i<count;i++)
		delete [] Ysend.data[i];
		delete [] Ysend.data;
		delete [] Y;
		delete [] X;
		delete [] phm_neighbours_index;


	}		//end of while


	delete [] resample_neighbours;


}		//end of if ( root != rank )


MPI_Barrier(MPI_COMM_WORLD);
cout << " I AM AT FINALIZE  rank = " << rank << "\n";
MPI_Finalize();

return 0;
}




void BuildX(double * X, struct Data4 * proximity, vector<struct Data3> & uniq_coord,struct Data5 & I,int i,struct Data2 * q,struct Data1 *p,int sample)
//void BuildX(double ** X, int iproc, struct Data4 * proximity, vector<struct Data3> & uniq_coord,struct Data5 & I,int i,struct Data2 * q,struct Data1 *p,int sample)
{
//int count;
int k,l,j;

	//This is to ensure that the first row of data is that of the point around which the emulator is centered
	X[0] = uniq_coord[proximity[i].uniq_id].x;
	X[1] = uniq_coord[proximity[i].uniq_id].y;

	for(k=0;k<Ndim;k++)
	X[k+2] = I.data[sample][k];

	X[Ndim+2] = p->spatial[i].H;
	l = 1;

	for(j=0;j< q->neighbour_count; j++)					//This loops over the neighbours of the centre to create data structure for building R
	{

		if ( proximity[i].uniq_id == proximity[i].uniq_id_neighbour[j] && sample == proximity[i].macro_neighbour[j] )
		continue;

		X[l*(Ndim+3)] = uniq_coord[proximity[i].uniq_id_neighbour[j]].x;
		X[l*(Ndim+3)+1] = uniq_coord[proximity[i].uniq_id_neighbour[j]].y;

		for(k=0;k<Ndim;k++)
		X[l*(Ndim+3)+k+2] = I.data[proximity[i].macro_neighbour[j]][k];

		X[l*(Ndim+3)+Ndim+2] =  p->spatial[proximity[i].spatial_neighbour[j]].H;

		l++;

	}


	if ( l > j )
	q->neighbour_count = q->neighbour_count+1;

}




//stdout the mean and variance information
void Send_To_Netezza(list<struct Data6> & Pipeout, int rank, int group, int fid)
{
time_t start,stop;
int count=0;

//int fid;
int j=0,i;
char pipename[25],output[300];
//Loading_file[25],stat_file[25];
list<struct Data6> :: iterator k;

/*
if ( *stat_count == 0 )
{
	sprintf(Loading_file,"./nzload.sh %d %d",rank,group);
	sprintf(pipename,"Netpipe%d",rank);
	system(Loading_file);
	fid = open(pipename,O_WRONLY);
}
*/

k = Pipeout.begin();

time(&start);

/*
while( k != Pipeout.end() )
{

	for(i=0;i < k->Return_count;i=i+4)
	{

		if ( k->mean_and_variance[i] != -1 )
		{
			sprintf(output,"%d %d %d %d %lf %lf\n",k->Sample+1,k->Uniq_record+1,int(k->mean_and_variance[i]),int(k->mean_and_variance[i+1]),k->mean_and_variance[i+2],k->mean_and_variance[i+3]);
			write(fid,output,strlen(output));
			fflush(stdout);
		}

		else
		continue;

	}
	delete [] k->mean_and_variance;
	k->mean_and_variance = NULL;
	k++;
	j++;
}
*/


/////////////////////////////////////// ADDED ON SEPTEMBER 1, 2012 /////////////////////////////////////////////////

char * buffer;
buffer = new char[1024];
buffer[0] = '\0';

while( k != Pipeout.end() )
{

	for(i=0;i < k->Return_count;i=i+4)
	{

		if ( strlen(output) > 100 )
		{
//			cout << " output length = " << strlen(output) << "\n";
			continue;
		}

		if ( k->mean_and_variance[i] != -1 )
		{

			sprintf(output,"%d %d %d %d %lf %lf\n",k->Sample+1,k->Uniq_record+1,int(k->mean_and_variance[i]),int(k->mean_and_variance[i+1]),k->mean_and_variance[i+2],k->mean_and_variance[i+3]);
			strcat(buffer,output);

//			if ( strlen(output) > 100 )
//			cout << " output length = " << strlen(output) << "\n";

			count++;

			if ( strlen(buffer) > 900 )
			{
				write(fid,buffer,strlen(buffer));
//				fflush(stdout);
				buffer[0] = '\0';
			}
		}

		else
		continue;

	}

	delete [] k->mean_and_variance;
	k->mean_and_variance = NULL;
	k++;
	j++;
}

if ( strlen(buffer) != 0 )
{
	write(fid,buffer,strlen(buffer));
//	fflush(stdout);

}
delete [] buffer;


////////////////////////////////////////////////////////     END      //////////////////////////////////////////////////////////

time(&stop);
cout << " TIME TAKEN = " << difftime(stop,start) << "\n";
}
