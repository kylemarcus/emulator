#ifndef __BUILD_EMULATOR__H
#define __BUILD_EMULATOR__H
#include <iostream>
#include <stdio.h>
#include <cmath>
#include "matrixop.h"
#include <list>
#include <map>
#include <utility>
#include <sstream>
#include <fstream>

using namespace std;

#define DEBUG(arg) cerr << " I AM HERE " << __LINE__ << " " <<  #arg << endl;

struct Data2
{
int neighbour_count;
double H;
};

struct Data1
{
int spatial_count;
struct Data2 * spatial;
};

struct Data3
{
double x,y;
};

struct Data4
{
int uniq_id;
int * spatial_neighbour, * uniq_id_neighbour, * macro_neighbour;
};

struct Data5
{
double **data;
};

//phm_count is the number of uniq_id neighbors a phm_id has.

class CombinePHM {

public:
static vector<int> phm_count_list;


private:
int phm_id,phm_count;
double total_weight;
double weighted_mean;

public:
CombinePHM() {}
explicit CombinePHM (int phm_id) : phm_id(phm_id), phm_count(0), total_weight(0), weighted_mean(0) {}
void update_mean(double phm_xcoord, double phm_ycoord, double uniq_xcoord, double uniq_ycoord, double mean);
int Count() { return phm_count;}
int get_phm_id() { return phm_id; }
void IncreaseCount() { phm_count++; }
double get_mean() {
if ( total_weight != 0 ) return weighted_mean/total_weight;
else return -1; }

};


struct DataList {
int resample;
list <CombinePHM> output;

};

void CombinePHM :: update_mean(double phm_xcoord, double phm_ycoord, double uniq_xcoord, double uniq_ycoord, double mean) {

  double dist = sqrt ( pow ( (phm_xcoord-uniq_xcoord), 2 ) + pow ( (phm_xcoord-uniq_xcoord), 2 ) );
  total_weight += dist;
  weighted_mean += mean*dist;
  phm_count++;
}

typedef pair<double,double> double_pair;
typedef pair<int,int> int_pair;

double newton_method(matrixop & R,matrixop & Rinv, matrixop & G, matrixop & beta, matrixop & rhs,int inewt, int Ndim, int Npts, struct Data5 & X,double * y, double * log_corr_len_guess, double * log_max_corr_len, double * log_min_corr_len);
void Mean_and_Variance(matrixop & G, matrixop & Rinveps, matrixop & R, matrixop & Rinv, matrixop & beta, struct Data5 & xsx, struct Data5 & ysy, int Nx, int Ny, int Ndim, double sigma, double * log_corr_len_guess, double ymax, int sample, string uniq_coord, int resample, vector<int> phm,DataList &  my_list, vector<double*> & phm_coord, map<string,double_pair> & x_y_coord, vector<int> & phm_count_list, ofstream & out_file);


double newton_method(matrixop & R,matrixop & Rinv, matrixop & G, matrixop & beta, matrixop & rhs, int inewt, int Ndim, int Npts, struct Data5 & X, double * y, double * log_corr_len_guess, double * log_max_corr_len, double * log_min_corr_len)
{
int i,j,k;
double temp = 0;
double sum=0;

double min_allowed_rcond;
min_allowed_rcond = pow(2,-40);

double * theta_guess = new double [Ndim];						//delete

R = matrixop(Npts,Npts);
Rinv = matrixop(Npts,Npts);

for(j=0;j<Ndim;j++)
theta_guess[j] = 0.5*exp(-2*log_corr_len_guess[j]);


/*cout <<"\nrank = " << rank << "  inewt = " << inewt << "   THETA GUESS     ::";
for(j=0;j<Ndim;j++)
cout << theta_guess[j] << " ";
cout << "\n";
*/

// Creating Y vector
matrixop Y(Npts,1);
for(i=0;i<Npts;i++)
Y.insert(y[i]);


if ( inewt == 0 )
{
G = matrixop(Npts,Ndim+1);

// G matrix
for(i=0;i<Npts;i++)
G.insert(1);

for(j=0;j<Ndim;j++)
{
	for(i=0;i<Npts;i++)
	G.insert(X.data[i][j]-X.data[0][j]);
}

}

/*
cout << "X.data =" << "\n";
for (i=0;i<Npts;i++) {
 for(j=0;j<Ndim;j++)
  cout << X.data[i][j] << " ";
 cout << "\n";
}
*/

// R matrix (Co-variance matrix)
for (i=0;i<Npts;i++)
{
	for(j=0;j<Npts;j++)
	{
		sum = 0;

		for(k=0;k<Ndim;k++)
		sum = sum + theta_guess[k] * -(pow( ( X.data[i][k] - X.data[j][k] ), 2 ));

		R.insert(exp(sum));
		Rinv.insert(exp(sum));							//Rinv because it will be inverted later

	}
}
//R.disp();

///////////////////////////////////////////////////////////////////    MODIFIED ON JULY 19   ///////////////////////////////////////////
double Rcond;
//double worstcaseleigmax,nugget;
Rcond = R.rcond();
/*
if ( Rcond < min_allowed_rcond )
{
	cout << "\n CHECK : ";
        worstcaseleigmax=Npts/(1+(Npts-1)*Rcond);
        nugget = worstcaseleigmax*(min_allowed_rcond-Rcond)/(1-min_allowed_rcond);

        for(j=0;j<Npts*Npts;j=j+Npts+1)
        R.change(j,( nugget + R.get(j) ) );
}
else
nugget = 0;

//Rinv = inv(R)
Rinv.inv();
*/

if ( Rcond < min_allowed_rcond )
Rinv.pinv();
else
Rinv.inv();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

matrixop Gtran_Rinv,Gtran_Rinv_G_inv,epsilon,mat;

//Gtran_Rinv = G'*Rinv
Gtran_Rinv.mul(G,'T',Rinv,'N',1);

//Gtran_Rinv_G_inv = inv(G'*Rinv*G)
Gtran_Rinv_G_inv.mul(Gtran_Rinv,'N',G,'N',1);
Gtran_Rinv_G_inv.inv();

mat.mul(Gtran_Rinv,'N',Y,'N',1);
beta.mul(Gtran_Rinv_G_inv,'N',mat,'N',1);
mat.mul(G,'N',beta,'N',1);

epsilon = matrixop( Y.r(),Y.c() );

for(i=0;i<(Y.r())*(Y.c());i++)
epsilon.insert( Y.get(i)-mat.get(i) );

rhs = epsilon;
rhs.div(R);		//This is same as R\epsilon or inverse(R)*epsilon

double unadj_var=0,Obj,loglike;

for(i=0;i<rhs.r()*rhs.c();i++)
unadj_var = unadj_var + rhs.get(i)*epsilon.get(i);

unadj_var = unadj_var/Npts;

//NEED A ROUTINE TO EVALUATE DETERMINANT
loglike = -0.5*(log(unadj_var)+log(abs(R.det()))/Npts);
Obj = -loglike;

if (inewt == 10)
return unadj_var;						//Check here what to return

double * dunadj_var = new double[Ndim];
double trace_Rinv_dR;
double *dObj = new double[Ndim];
int idim;

matrixop dRdi,dR_rhs,dbeta,depsilon;

matrixop mat1;
matrixop mat2;

dR_rhs = matrixop(Npts,Ndim);
depsilon = matrixop(G.r(),Ndim);

for(idim=0;idim<Ndim;idim++)
{
	dRdi = matrixop(Npts,Npts);

	for(j=0;j<Npts;j++)
	{
		for(k=0;k<Npts;k++)
		{
			temp = -1 * R.get(j*Npts+k) * pow( (X.data[j][idim] - X.data[k][idim]),2);
			dRdi.insert(temp);
		}
	}

	mat2.mul(dRdi,'N',rhs,'N',1);

	for(i=0;i<mat2.r()*mat2.c();i++)
	dR_rhs.insert(mat2.get(i));

	mat1.mul(Gtran_Rinv,'N',mat2,'N',1);
	dbeta.mul(Gtran_Rinv_G_inv,'N',mat1,'N',-1);
	mat.mul(G,'N',dbeta,'N',-1);

	for(i=0;i<mat.r()*mat.c();i++)
	depsilon.insert(mat.get(i));

	mat1.mul(rhs,'T',mat,'N',2);
	mat.mul(rhs,'T',mat2,'N',1);
	dunadj_var[idim] = ( mat1.get(0) - mat.get(0) )/Npts;

	temp = 0;
	for(j=0;j<Rinv.r()*Rinv.c();j++)
	temp = temp + Rinv.get(j)*dRdi.get(j);
	trace_Rinv_dR = temp;

	dObj[idim] = 0.5*(dunadj_var[idim]/unadj_var+trace_Rinv_dR/Npts);

}

matrixop dR_rhs_minus_depsilon(depsilon.r(),depsilon.c());

for(i=0;i<depsilon.r()*depsilon.c();i++)
dR_rhs_minus_depsilon.insert(dR_rhs.get(i)-depsilon.get(i));

matrixop Rinv_dR_rhs_minus_depsilon;
Rinv_dR_rhs_minus_depsilon.mul(Rinv,'N',dR_rhs_minus_depsilon,'N',1);






dRdi = matrixop(Npts,Npts);
idim = 0;

for(j=0;j<Npts;j++)
{
	for(k=0;k<Npts;k++)
        dRdi.insert(R.get(j*Npts+k)* -pow( (X.data[j][idim] - X.data[k][idim]),2) );
}

matrixop Rinv_dRdi;
Rinv_dRdi.mul(Rinv,'N',dRdi,'N',1);

matrixop temp4(dR_rhs_minus_depsilon.r(),idim+1);

for(i=0;i<dR_rhs_minus_depsilon.r();i++)
temp4.insert(dR_rhs_minus_depsilon.get(i));

double d2unadj_vardidj = 0;
matrixop temp2,d2Rdidj,temp5,temp1,d2beta,d2epsilon,d2Obj;
matrixop dRdj,temp3,Rinv_dRdj;

d2Obj = matrixop(Ndim,Ndim);
//temp1 is same as temp in Keith's code
int jdim;

for(idim=0;idim<Ndim;idim++)
{
	temp2 = matrixop(Rinv_dR_rhs_minus_depsilon.r(),1);
	d2Rdidj = matrixop(Npts,Npts);

	for(i=0;i<Rinv_dR_rhs_minus_depsilon.r();i++)
	temp2.insert(Rinv_dR_rhs_minus_depsilon.get( idim*( Rinv_dR_rhs_minus_depsilon.r() )+i ));

	jdim = idim;

	for(j=0;j<Npts;j++)
	{
       		for(k=0;k<Npts;k++)
	        d2Rdidj.insert(dRdi.get(j*Npts+k)* -pow( (X.data[j][jdim] - X.data[k][jdim]),2) );
	}

	temp5.mul(d2Rdidj,'N',rhs,'N',-1);
	mat.mul(rhs,'T',temp5,'N',1);
	d2unadj_vardidj = mat.get(0);
	mat.mul(dRdi,'N',temp2,'N',2);

	for(i=0;i<temp5.r()*temp5.c();i++)
	temp5.change(i,temp5.get(i)+mat.get(i));

	temp1.mul(Gtran_Rinv,'N',temp5,'N',1);
	d2beta.mul(Gtran_Rinv_G_inv,'N',temp1,'N',1);
	d2epsilon.mul(G,'N',d2beta,'N',-1);
	mat1.mul(rhs,'T',d2epsilon,'N',2);
	mat2.mul(temp2,'T',temp4,'N',2);
	d2unadj_vardidj = ( d2unadj_vardidj + mat1.get(0) + mat2.get(0) )/Npts;

	temp = 0;
	k = 0;

	for(i=0;i<Rinv_dRdi.r()*Rinv_dRdi.c();i++)
	{
		j = i%Rinv_dRdi.r();
		k = int(i/Rinv_dRdi.r());
		temp = temp + Rinv_dRdi.get(i)*Rinv_dRdi.get(j*Rinv_dRdi.r()+k);
	}

	temp = temp*-1;
	for(i=0;i<Rinv.r()*Rinv.c();i++)
	temp = temp + Rinv.get(i)*d2Rdidj.get(i);
	temp = temp/Npts;
	temp = 0.5 * (-1*pow((dunadj_var[idim]/unadj_var),2) + d2unadj_vardidj/unadj_var + temp);
	d2Obj.change((idim*Ndim+idim),temp);

	temp3 = matrixop(Rinv_dR_rhs_minus_depsilon.r(),1);

	for(jdim=Ndim-1;jdim>idim;jdim--)
	{

		dRdj = matrixop(Npts,Npts);
		d2Rdidj = matrixop(Npts,Npts);

		for(j=0;j<Npts;j++)
		{
	       		for(k=0;k<Npts;k++)
		        dRdj.insert(R.get(j*Npts+k)* -pow( (X.data[j][jdim] - X.data[k][jdim]),2) );
		}

		for(j=0;j<Npts;j++)
		{
	       		for(k=0;k<Npts;k++)
		        d2Rdidj.insert(dRdi.get(j*Npts+k)* -pow( (X.data[j][jdim] - X.data[k][jdim]),2) );
		}


		for(i=0;i<Rinv_dR_rhs_minus_depsilon.r();i++)
		temp3.change(i,Rinv_dR_rhs_minus_depsilon.get(jdim*Rinv_dR_rhs_minus_depsilon.r()+i));

		Rinv_dRdj.mul(Rinv,'N',dRdj,'N',1);
		temp5.mul(d2Rdidj,'N',rhs,'N',-1);
		mat.mul(rhs,'T',temp5,'N',1);
		d2unadj_vardidj = mat.get(0);
		mat.mul(dRdi,'N',temp3,'N',1);
		mat1.mul(dRdj,'N',temp2,'N',1);
		for(i=0;i<temp5.r()*temp5.c();i++)
		temp5.change(i,(temp5.get(i)+mat.get(i)+mat1.get(i)) );
		temp1.mul(Gtran_Rinv,'N',temp5,'N',1);					//temp1 is same as temp in Keith's code
		d2beta.mul(Gtran_Rinv_G_inv,'N',temp1,'N',1);
		d2epsilon.mul(G,'N',d2beta,'N',-1);
		for(i=0;i<dR_rhs_minus_depsilon.r();i++)
		temp4.change(i,dR_rhs_minus_depsilon.get( (jdim*dR_rhs_minus_depsilon.r()) + i) );
		mat.mul(rhs,'T',d2epsilon,'N',1);
		temp = mat.get(0);
//		cout << temp << "\n ";
		mat.mul(temp2,'T',temp4,'N',1);
		d2unadj_vardidj = (d2unadj_vardidj + 2*( temp + mat.get(0) ) )/Npts;
//		cout << "d2unadj_vardidj = " << d2unadj_vardidj << "\n";
		temp = 0;
	        k = 0;

		for(i=0;i<Rinv_dRdi.r()*Rinv_dRdi.c();i++)
		{
			j = i%Rinv_dRdj.r();
			k = int(i/Rinv_dRdj.r());
			temp = temp + Rinv_dRdj.get(i)*Rinv_dRdi.get(j*Rinv_dRdj.r()+k);
		}

	        temp = temp*-1;
	        for(i=0;i<Rinv.r()*Rinv.c();i++)
	        temp = temp + Rinv.get(i)*d2Rdidj.get(i);
	        temp = d2unadj_vardidj/unadj_var + temp/Npts;
	        temp = 0.5*(-1*(dunadj_var[idim]/unadj_var * dunadj_var[jdim]/unadj_var) + temp);
	        d2Obj.change((idim*Ndim+jdim),temp);
	        d2Obj.change((jdim*Ndim+idim),temp);

	}


	if(Ndim>idim+1)
	{
	        dRdi=dRdj;
	        Rinv_dRdi=Rinv_dRdj;
	}

}


//dObj.disp();

matrixop dtheta_to_dlog_corr_len(Ndim,Ndim),dir;
j = 0;
k = 0;
for(i=0;i<Ndim*Ndim;i++)
{
	if ( i == k )
	{
		dtheta_to_dlog_corr_len.insert(-2*theta_guess[j++]);
		k = k+Ndim+1;
	}

	else
	dtheta_to_dlog_corr_len.insert(0);
}

mat2 = matrixop(Ndim,1);
for(i=0;i<Ndim;i++)
mat2.insert(dObj[i]);

mat.mul(dtheta_to_dlog_corr_len,'N',mat2,'N',1);
//dObj = mat;
mat2 = mat;
mat.mul(dtheta_to_dlog_corr_len,'N',d2Obj,'N',1);   		//mat2 is dObj from here
mat1.mul(mat,'N',dtheta_to_dlog_corr_len,'N',1);
d2Obj = mat1;

dir = mat2;

dir.div(d2Obj);
mat.mul(mat2,'T',dir,'N',1);

temp = mat.get(0);
temp = temp/abs(temp);

for(i=0;i<dir.r()*dir.c();i++)
dir.change(i,dir.get(i)*-1*temp);

double * stopfact = new double [dir.r()*dir.c()];
double * ifdirzero = new double [dir.r()*dir.c()];

int finite[Ndim];
k = 0;
for(i=0;i<dir.r()*dir.c();i++)
{
	ifdirzero[i] = (dir.get(i) == 0) ? 1:0;

	if ( dir.get(i)+ifdirzero[i] != 0 )
	{
		stopfact[i] = ( ( (log_max_corr_len[i]-log_corr_len_guess[i])*( (dir.get(i)>0) ? 1:0 ) + (log_corr_len_guess[i]-log_min_corr_len[i])*( (dir.get(i)>0) ? 0:1 ) )/abs( dir.get(i)+ifdirzero[i] ) )*(1-ifdirzero[i]);
		finite[i] = 1;
	}

	if ( dir.get(i)+ifdirzero[i] == 0 )
	{
		stopfact[i] = 10^12;
		finite[i] = 0;
	}

	if ( stopfact[i] == 0 )
	k= 1;

	if ( stopfact[i] < 0 )
	finite[i] = 0;
}


if ( k== 1)
{
	for(i=0;i<dir.r()*dir.c();i++)
	dir.change(i,dir.get(i)*0.05);
}


for(i=0;i<Ndim;i++)
dir.change(i,dir.get(i)*finite[i]);

double stopfact2;

/*
for(i=0;i<Ndim;i++)
{
k = (dir.get(i) == 0) ? 1:0;

stopfact[i] = stopfact[i]*(1-k)+k;
}
*/

for(i=0;i<Ndim;i++)
{
if ( i == 0 )
stopfact2 = stopfact[i];

else
stopfact2 = ( stopfact2 < stopfact[i] ) ? stopfact2:stopfact[i];

}

stopfact2 = (stopfact2>1) ? 1:stopfact2;

for(i=0;i<Ndim;i++)
dir.change(i,dir.get(i)*stopfact2);

//dir.disp();
//cout <<" stopfact2 = " << stopfact2<< "\n";

for(i=0;i<Ndim;i++)
log_corr_len_guess[i] = log_corr_len_guess[i] + dir.get(i);

delete [] theta_guess;					/////////////////////////////////////// ADDED ON JULY 18     4:46pm  ////////////////////////////
return unadj_var;


}

void Mean_and_Variance(matrixop & G, matrixop & Rinveps, matrixop & R, matrixop & Rinv, matrixop & beta, struct Data5 & xsx, struct Data5 & ysy, int Nx, int Ny, int Ndim, double sigma, double * log_corr_len_guess, double ymax, int sample, string uniq_coord, int resample, vector<int> phm, DataList &  my_list, vector<double*> & phm_coord, map<string,double_pair> & x_y_coord, vector<int> & phm_count_list, ofstream & out_file)

{
int i,j,k;
matrixop g(Nx,Ndim+1),rTran(Nx,Ny),rTran_Rinv,temp;
matrixop mat,mat1,mat2,mat3;
double theta_guess[Ndim];
double sum;

for(j=0;j<Ndim;j++)
theta_guess[j] = 0.5*exp(-2*log_corr_len_guess[j]);

for (i=0;i<Ny;i++)
{
        for(j=0;j<Nx;j++)
        {
                sum = 0;

                for(k=0;k<Ndim;k++)
                sum = sum + theta_guess[k] * -(pow( ( xsx.data[j][k] - ysy.data[i][k] ), 2 ));

                rTran.insert(exp(sum));

        }
}


for(i=0;i<Nx;i++)
g.insert(1);

for(j=0;j<Ndim;j++)
{
        for(i=0;i<Nx;i++)
        g.insert(xsx.data[i][j]-ysy.data[0][j]);
}


rTran_Rinv.mul(rTran,'N',Rinv,'N',1);
//rTran_Rinv.div(R);


mat.mul(rTran_Rinv,'N',G,'N',1);
temp = matrixop(g.r(),g.c());

for(i=0;i<g.r()*g.c();i++)
temp.insert(g.get(i)-mat.get(i));

mat.mul(g,'N',beta,'N',1);

mat1.mul(rTran,'N',Rinveps,'N',1);

//Variance calculation
mat2 = G;
mat2.mul(Rinv,'N',G,'N',1);
mat3.mul(G,'T',mat2,'N',1);
mat3.inv();
mat2.mul(temp,'N',mat3,'N',1);

//mat2.mul(rTran_Rinv,'N',rTran,'N',1);
double * sum1,* sum2;
int num;

sum1 = new double [temp.r()];

for(i=0;i<temp.r();i++)
{
	sum1[i] = 0;
	for(j=0;j<temp.c();j++)
	{
		num = j*temp.r();
		sum1[i] = ( mat2.get(num+i)*temp.get(num+i) ) + sum1[i];

	}
}


mat3.mul(mat2,'N',temp,'T',1);
sum2 = new double [rTran.r()];

for(i=0;i<rTran.r();i++)
{
	sum2[i] = 0;

	for(j=0;j< rTran.c();j++)
	{
		num = j*rTran.r();
		sum2[i] = rTran_Rinv.get(num+i)*rTran.get(num+i) + sum2[i];

	}
}



//mat2.mul(rTran_Rinv,'N',rTran,'T',1);

//list<DataStruct> :: iterator phm_iter;

list<CombinePHM> :: iterator phm_iter;
char str[100];
stringstream ss;

for(i=0;i<mat.r();i++)
{

sum  = mat.get(i) + mat1.get(i);

//if ( sum <= 1.5*ymax )  {
   double mean = ( sum > 0 ) ? sum:0;
   mean = (mean<=ymax*1.5) ? mean : -1;
   int phm_id = phm[i];

   bool flag = 0;
   phm_iter = my_list.output.begin();

   while ( phm_id >= phm_iter->get_phm_id() && phm_iter != my_list.output.end() ) {

	if ( phm_id == phm_iter->get_phm_id() ) {
		flag = 1;
		break;
	}//end of if
	phm_iter++;

    }//end of while

   if ( !flag ) phm_iter = my_list.output.insert( phm_iter, CombinePHM(phm_id) );

   if ( mean >= 0 )
     phm_iter->update_mean( phm_coord[phm_id][0], phm_coord[phm_id][1], x_y_coord[uniq_coord].first,  x_y_coord[uniq_coord].second, mean);
   else  phm_iter->IncreaseCount();

   if ( phm_iter->Count() == phm_count_list[phm_id] ) {
	sprintf(str,"%d_%d",resample,phm_id);
//	cout << str << " " << sample << " " << phm_iter->get_mean() << endl;
	ss << str << " " << sample << " " << phm_iter->get_mean() << "\n";
//      cout << sample << " " << resample << " " << phm_id << " " << phm_iter->get_mean() << endl;
      phm_iter = my_list.output.erase(phm_iter);

   }//end of if

   if ( ss.str().length() > 1000 ) {
//	cout << ss.str().c_str();
       out_file.write(ss.str().c_str(),ss.str().length());
       ss.str("");
   }


}

if ( ss.str().length() > 0 )
out_file.write(ss.str().c_str(),ss.str().length());
//out_file << ss.str().c_str();

//}

delete [] sum1;
delete [] sum2;

}


#endif
