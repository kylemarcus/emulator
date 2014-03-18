#include "mpi.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <fstream>

using namespace std;

int main( int argc, char * argv[] )
{

int rank,num_procs;
MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD,&rank);
MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
int sample_count;

if (! argv[1] ){
cout << " Please provide the number of pile height records that need to be downsampled" << endl;
}

sample_count  = atoi(argv[1]);
cout << "sample_count = " << sample_count << endl;
//char file_name[50];

//cout << file_name << endl;

if ( rank == 0 )
{

 int recv_flag[num_procs-1];
 int recv_stat[num_procs-1];
 MPI_Request req[num_procs-1],req_recv[num_procs-1];
 MPI_Status stat[num_procs-1];
 int count = 0;

 for(int it=1; it<=num_procs; it++)
 recv_stat[it-1] = 0;

 for(int it=1; it<num_procs; it++)
 MPI_Irecv(&recv_stat[it-1],1,MPI_INT,it,it+num_procs,MPI_COMM_WORLD,&req_recv[it-1]);

 for(int it=1; it<num_procs; it++)
 MPI_Send(&it,1,MPI_INT,it,it,MPI_COMM_WORLD);
 count = num_procs - 1;

 while ( count < sample_count )
 {
   for(int iproc=1;iproc<num_procs;iproc++)
    {
	MPI_Test(&req_recv[iproc-1],&recv_flag[iproc-1],&stat[iproc-1]);
	if ( recv_flag[iproc-1] && count < sample_count )
	    {
		recv_flag[iproc-1] = 0;
		count++;
                MPI_Irecv(&recv_stat[iproc-1],1,MPI_INT,iproc,iproc+num_procs,MPI_COMM_WORLD,&req_recv[iproc-1]);
		MPI_Send(&count,1,MPI_INT,iproc,iproc,MPI_COMM_WORLD);
	    }
     }
 }

 int final_flag=0;
 for(int it=1; it<num_procs; it++)
 MPI_Send(&final_flag,1,MPI_INT,it,it,MPI_COMM_WORLD);


}


if ( rank != 0 )
{
 string file;
 vector<string> file_list;
// ifstream my_file("file_list_4");
 ifstream my_file(argv[2]);
 if (my_file.good()) cout << " opened file " << endl;
 else cout << " file does not exist " << endl;

 while ( my_file >> file ) file_list.push_back(file);
 my_file.close();
 sample_count = file_list.size();


MPI_Status stat;
int sample_num=0;
char filename[100];
MPI_Recv(&sample_num,1,MPI_INT,0,rank,MPI_COMM_WORLD,&stat);

while ( sample_num ) {


       	MPI_Send(&sample_num,1,MPI_INT,0,num_procs+rank,MPI_COMM_WORLD);
//	sprintf(filename,"./down_sample.sh %d",sample_num);
//	sprintf(filename,"./phm_count.sh %d",sample_num);
//	sprintf(filename,"./map3.sh %d",sample_num);
	cout << "   filename = " << file_list[sample_num-1].c_str() << endl;
	sprintf(filename,"./map6.sh /panasas/scratch/shivaswa/output/%s",file_list[sample_num-1].c_str());
       	cout << " rank = " << rank << "    sample = " << sample_num << "  filename = " << filename <<endl;
	system(filename);
	MPI_Recv(&sample_num,1,MPI_INT,0,rank,MPI_COMM_WORLD,&stat);

  }
}

MPI_Barrier(MPI_COMM_WORLD);

MPI_Finalize();
return 0;
}



