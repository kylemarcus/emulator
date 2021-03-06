#include <iostream>
#include <vector>
#include "matrixop.h"
#include "build_emulator.h"
#include <fstream>
#include <list>
#include <map>
#include <sstream>
#include <stdio.h>
#include <time.h>

using namespace std;

const int RADIUS = 100;
const int MAX_SIZE = 200;
const double MACRO_SCALED_RADIUS = 0.2;
const int Ndim = 4;
#define DEBUG(arg) cerr << " I AM HERE " << __LINE__ << " " <<  #arg << endl;

void display(vector<double*> & my_data);
void parse ( vector<double> & param, double **X, double *Y);
void corr(double *log_corr_len_guess, double *log_min_corr_len, double *log_max_corr_len, double *lmin );


int main() {
time_t start,stop;
string input_stream;
vector<string> input;
input.reserve(MAX_SIZE);
//double Y[MAX_SIZE];
Data5 X;
double * Y;
double max,min,log_min_corr_len[Ndim+2],log_max_corr_len[Ndim+2], lmin[Ndim],log_corr_len_guess[Ndim+2], sigma;

//Read uniq_coords.txt
FILE *fp;
stringstream ss;
typedef pair<double,double> double_pair;
typedef pair<int,int> int_pair;


for(int i = 0; i<Ndim;i++)
{
    getline(cin,input_stream);
    sscanf(input_stream.c_str(),"%lf %lf",&max,&min);
    lmin[i] =  0.25*MACRO_SCALED_RADIUS*(max-min);
}


//read Nx_Ny.txt
double Xstart,Xend,Ystart,Yend;
vector<int_pair> Nx_Ny;
fp = fopen("/user/shivaswa/my_hadoop/Nx_Ny_group.txt","r");
//fgets(store,100,fp);
fscanf(fp,"%lf %lf %lf %lf", &Xstart,&Xend,&Ystart,&Yend);
while ( !feof(fp) ) {
int Nx, Ny, group;
fscanf(fp,"%d %d %d",&group,&Nx,&Ny);
int_pair temp = make_pair(Nx,Ny);
Nx_Ny.push_back(temp);
}
fclose(fp);


fp = fopen("/user/shivaswa/my_hadoop/uniq_coords.txt","r");
map<string,double_pair> x_y_coord;
while ( !feof(fp) )
{
int row,col,group,Nx,Ny;
double x,y;
fscanf(fp,"%d_%d_%d",&group,&row,&col);
Nx = Nx_Ny[group].first;
Ny = Nx_Ny[group].second;
ss << group << "_" << row << "_" << col;
row = row - 1;
col = col - 1;
x = (2*col + 0.5)*( Xend - Xstart )/(2*Nx) + Xstart;
y = (2*row + 0.5)*( Yend - Ystart )/(2*Ny) + Ystart;
double_pair x_y_pair = make_pair(x,y);
x_y_coord[ss.str()] = x_y_pair;
//cout << ss.str() << "  " << x_y_coord[ss.str()].first << "  " << x_y_coord[ss.str()].second << " " << group << " " << row << " " << col << endl;
ss.str("");

}

fclose(fp);


// Read resamples file
double *read_file;
vector<double*> resamples;
int garb;
fp = fopen("/user/shivaswa/my_hadoop/resamples.txt","r");
do {
read_file = new double[Ndim+1];
fscanf(fp,"%d,%lf,%lf,%lf,%lf,%lf",&garb,&read_file[0],&read_file[1],&read_file[2],&read_file[3],&read_file[4]);
resamples.push_back(read_file);
} while ( !feof(fp) );
fclose(fp);

vector<double*> phm;
fp = fopen("/user/shivaswa/my_hadoop/montserrat_take2_vol_dir_bed_int.phm","r");
do {
read_file = new double[2];
fscanf(fp,"%d %lf %lf",&garb,&read_file[0],&read_file[1]);
phm.push_back(read_file);
} while ( !feof(fp) );
fclose(fp);

ofstream myfile;
myfile.open("/user/shivaswa/my_hadoop/BLAH.txt",ofstream::out);
myfile.close();

resamples.pop_back();
//display(resamples);
vector<int> res_neigh;
int read_int;
double read_double;

int sample,count=0;
//bool next_sample_flag=0;
vector<int> phm_neigh;
DataStruct output;
list< DataList > my_list;
list< DataList > :: iterator my_iter;
list< DataStruct > :: iterator output_iter, second_iter;



do {

    getline(cin,input_stream);
//    cout << input_stream << endl;

    if ( !input_stream.compare("SAMPLE") ) {
	cin.clear();
	cin >> sample;
	cin >> input_stream;
	cin >> input_stream;
//	cout << " sample = " << sample << endl;
	myfile.open("/user/shivaswa/my_hadoop/BLAH.txt",ofstream::app);
	myfile << "sample = " << sample << endl;
	myfile.close();
//	DEBUG()
//	cout << " NEXT SAMPLE = " << sample << "<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>" << endl;

	if ( my_list.size() ) {
	    double x,y,dist;
	    string uniq_id;
	    vector<double> dist_vec;
	    for(my_iter=my_list.begin(); my_iter != my_list.end(); my_iter++ ) {
		output_iter = my_iter->output.begin();
		second_iter = output_iter;
		double dist_sum=0;
         	int phm_id = output_iter->phm_id;
		while ( output_iter != my_iter->output.end() ) {
		    uniq_id = output_iter->uniq_id;
		    x = x_y_coord[uniq_id].first;
		    y = x_y_coord[uniq_id].second;

		    if ( phm_id == output_iter->phm_id ) {
		      dist =  sqrt( pow( (phm[phm_id][0]-x),2 ) + pow( (phm[phm_id][1]-y),2 ) );
		      dist_sum = dist_sum + dist;
		      dist_vec.push_back(dist);
		    }
		    else {
		      double weighted_mean=0;
		      int j=0;
		      while ( second_iter != output_iter )   {
			 weighted_mean =  weighted_mean + second_iter->mean*dist_vec[j]/dist_sum;
			 second_iter++;
			 j++;
		      }

		      if (second_iter != output_iter) {
			 DEBUG("<<<<<<<<< TROUBLE  AT LINE  >>>>>>>>>>")
			 exit(1); }

		      stringstream ss_out;
		      ss_out << my_iter->resample << "_" << phm_id;
		      cout << ss_out.str() << " " << sample << " " << weighted_mean << endl;

		      phm_id = output_iter->phm_id;
		      dist_sum = sqrt( pow( (phm[phm_id][0]-x),2 ) + pow( (phm[phm_id][0]-y),2 ) );
		      dist_vec.clear();
		      dist_vec.push_back(dist_sum);
		    }

		    output_iter++;

		}//end of while
	   }//end of for
	}//end of if


	for(my_iter=my_list.begin(); my_iter != my_list.end(); my_iter++ ) my_iter->output.clear();
	my_list.clear();
       }

    if ( !input_stream.compare("RESAMPLES") ) {
	cin.clear();
//	cin >> read_int;
	while ( cin >> read_int ) {
//		cout << " READ = " << read_int << endl;
		res_neigh.push_back(read_int); }
		cin.clear();
	}


    if ( !input_stream.compare("NEXT") )
    {
	time(&start);
	string uniq_coord;
	cin >> uniq_coord;
//	cout << " uniq_coord = " << uniq_coord << endl;
	matrixop R,beta,rhs,G,Rinv;
	vector<double> param;
	while ( cin >> read_double ) param.push_back(read_double);
	cin.clear();

	int Npts = int(param.size())/(Ndim+3);
	Y = new double[Npts];
	X.data = new double*[Npts];
	for (int idata=0; idata< Npts; idata++)
	X.data[idata] = new double[Ndim+2];
	parse( param, X.data, Y );
//	cout << "Npts = " << Npts << endl;
//	Npts = parse( input, X.data, Y );
	time(&stop);

	cout << count++ << " " << " time to read = " << difftime(stop,start) << endl;
	corr(log_corr_len_guess, log_min_corr_len, log_max_corr_len, lmin );
	for (int inewt=0;inewt<10;inewt++)
	sigma = newton_method(R,Rinv,G,beta,rhs,inewt,Ndim+2,Npts,X,Y,log_corr_len_guess,log_max_corr_len,log_min_corr_len);

	cin >> input_stream;
	if ( !input_stream.compare("STOP") ) cin >> input_stream;

	 if ( ! input_stream.compare("PHM") ) {
	  phm_neigh.clear();
	  while ( cin >> read_int ) phm_neigh.push_back(read_int);
	  cin.clear();
	 }

	 if ( !phm_neigh.size() ) cerr << "PHM neighbors missing " << endl;

	double ymax;
	for(int i=0;i<Npts;i++)
        {
             if ( ymax < Y[i] )
             ymax = Y[i];
        }

	vector<int> :: iterator iter = res_neigh.begin();
	vector<Data5*> Res_X;

	int Mpts = phm_neigh.size();

	for(int ires=0;ires<int(res_neigh.size());ires++) {
	    double * mark_res = resamples[ires];
	    Data5 * X_r = new Data5;
	    X_r->data = new double*[Mpts];
	    for (int iall=0;iall<Mpts;iall++)
	    X_r->data[iall] = new double[Ndim+2];

	    for (int iphm=0;iphm<int(phm_neigh.size());iphm++) {
		double * mark_phm = phm[iphm];
		X_r->data[iphm][0] = *mark_phm++;
		X_r->data[iphm][1] = *mark_phm;
		for(int imark=0;imark<Ndim;imark++)
		X_r->data[iphm][imark+2] = mark_res[imark];
	    }//end of for

	   Res_X.push_back(X_r);
	}//end of for

	for (int iresample=0;iresample< int(res_neigh.size()); iresample++) {

	    int curr_resample = res_neigh[iresample];
	    DataList temp_datalist;

	    bool resample_flag=0;
	    my_iter = my_list.begin();

	    while ( my_iter != my_list.end() ) {
		if ( curr_resample > my_iter->resample ) my_iter++;
		else if ( curr_resample == my_iter->resample ) { resample_flag = 1;  break;}
		else break;
	    }

	    if ( !resample_flag ){
	        temp_datalist.resample = curr_resample;
		my_iter = my_list.insert(my_iter, temp_datalist );
	    }

//	  cout << " SENDING RESAMPLE = " << my_iter->resample << endl;
	  Mean_and_Variance(G, rhs, R, Rinv, beta, *Res_X[iresample], X, Mpts, Npts, Ndim+2, sigma, log_corr_len_guess, ymax, sample,uniq_coord,res_neigh[iresample], phm_neigh, *my_iter);
	  for(int i=0;i<Mpts;i++) delete [] Res_X[iresample]->data[i];
	  delete [] Res_X[iresample]->data;
	}


	Res_X.clear();

	for (int idata=0; idata< Npts; idata++)
	delete X.data[idata];
	delete [] X.data;
	delete [] Y;

    }


//	cout << sigma << endl;
} while ( input_stream.length() );

return 0;

}



void parse(vector<double> & param, double **X, double * Y )
{
 vector<double> :: iterator iter = param.begin();
 int Npts = int(param.size())/(Ndim+3);

 for(int i=0;i<Npts;i++) {
  for(int j=0;j<Ndim+2;j++) X[i][j] = *iter++;
  Y[i] = *iter++;
//  iter++;
//  cout << " Y = " << Y[i];
 }

}


void corr(double *log_corr_len_guess, double *log_min_corr_len, double *log_max_corr_len, double *lmin ) {

	log_min_corr_len[0] = log(RADIUS/4);
	log_min_corr_len[1] = log(RADIUS/4);
	log_max_corr_len[0] = log(2*RADIUS);
	log_max_corr_len[1] = log(2*RADIUS);

	//DEBUG(lmin)

	for(int i=0;i<Ndim;i++)
	{
	        log_min_corr_len[i+2] = log(lmin[i]);
	        log_max_corr_len[i+2] = log(8*lmin[i]);
	}

	log_corr_len_guess[0] = log(RADIUS);
	log_corr_len_guess[1] = log(RADIUS);

	for(int i=0;i<Ndim;i++)
	log_corr_len_guess[i+2] = log(4*lmin[i]);


}





void display( vector<double*> & my_data)
{
for (int i=0; i< int(my_data.size()); i++){
	for (int j=0;j < Ndim+1 ; j++)
	    cout << my_data[i][j] << " ";
	cout << endl;
}

}

