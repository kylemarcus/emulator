#include <iostream>
#include <vector>
#include "matrixop.h"
#include "build_emulator.h"
#include <fstream>
#include <list>

using namespace std;

const int RADIUS = 100;
const int MAX_SIZE = 200;
const double MACRO_SCALED_RADIUS = 0.2;
const int Ndim = 4;
#define DEBUG(arg) cout << " I AM HERE " << __LINE__ << " " <<  #arg << endl;

void display(vector<double*> & my_data);
void parse ( vector<double> & param, double **X, double *Y);
void corr(double *log_corr_len_guess, double *log_min_corr_len, double *log_max_corr_len, double *lmin );

struct DataStruct {
int uniq_id,phm_id;
double mean;

DataStruct() {}
DataStruct(int arg1, int arg2, double arg3) : uniq_id(arg1), phm_id(arg2), mean(arg3) {}
~DataStruct() {}
};



int main() {

string input_stream;
vector<string> input;
input.reserve(MAX_SIZE);
//double Y[MAX_SIZE];
Data5 X;
double * Y;
double max,min,log_min_corr_len[Ndim+2],log_max_corr_len[Ndim+2], lmin[Ndim],log_corr_len_guess[Ndim+2], sigma;


// Read resamples file
double *read_file;
vector<double*> resamples;
int garb;
FILE *fp;
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

int count=0;
ofstream myfile;
myfile.open("/user/shivaswa/my_hadoop/BLAH.txt",ofstream::out);
myfile.close();

resamples.pop_back();
//display(resamples);
vector<int> res_neigh;
int read_int;
double read_double;

for(int i = 0; i<Ndim;i++)
{
    getline(cin,input_stream);
    sscanf(input_stream.c_str(),"%lf %lf",&max,&min);
    lmin[i] =  0.25*MACRO_SCALED_RADIUS*(max-min);
}


int sample;
vector<int> phm_neigh;
DataStruct output;
list< list< DataStruct> > output_list;
list< list< DataStruct> > :: iterator output_iter;

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
	cout << count++ << " ";
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

	corr(log_corr_len_guess, log_min_corr_len, log_max_corr_len, lmin );
//	for (int inewt=0;inewt<10;inewt++)
//	sigma = newton_method(R,Rinv,G,beta,rhs,inewt,Ndim+2,Npts,X,Y,log_corr_len_guess,log_max_corr_len,log_min_corr_len);

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
//	while ( iter != res_neigh.end() ) {
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
	}//end of while

	for (int iresample=0;iresample< int(res_neigh.size()); iresample++) {
//	  Mean_and_Variance(G, rhs, R, Rinv, beta, *Res_X[iresample], X, Mpts, Npts, Ndim+2, sigma, log_corr_len_guess, ymax, sample,uniq_coord,res_neigh[iresample], phm_neigh);
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

