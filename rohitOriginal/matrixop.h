#ifndef _MATRIXOP_H
#define _MATRIXOP_H
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

extern "C"{

double dgemm_(char * transa, char * transb, int *M,int *N,int *K,double *alpha,double * A,int *lda,double * B,int *ldb,int *beta,double * C,int *ldc);
void dgetrf_(int * M, int * N, double * A, int * lda, int * INIV,int * INFO);
void dgecon_(char * NORM, int * M, double * A, int * lda, double * ANORM ,double * RCOND, double * WORK, int * IWORK, int * INFO);
void dgetri_(int * row, double * A,int * lda,int * IPIV, double *WORK, int *LWORK, int *INFO );
void dgesv_(int * row,int * NRHS,double * A,int * lda,int * IPIV,double * B,int * ldb,int * INFO );
void dgesvd_(char * JOBU, char * JOBVT, int * M, int * N, double * A, int * LDA, double * S, double * U, int * LDU, double * VT, int * LDVT, double * WORK, int * LWORK, int * INFO );

}

using namespace std;


class matrixop
{

private:
int row,col,k;
double *p;

public:
matrixop();
~matrixop();
matrixop(const matrixop & mat);
matrixop(int m, int n);
matrixop & operator=(const matrixop & mat);

int r();
int c();
void disp();
void insert(double x);
double get(int i);
void change(int i, double num);
void mul(matrixop & matA, char transA, matrixop & matB, char transB, double alpha);
double rcond();
void inv();
void div(matrixop & mat);            //vector.div(matrix) or b.div(A) in order to do A\b, note that b i.e. vector gets modified
double det();
void pinv();

//matrixop & operator*(const matrixop & mat);

};


matrixop :: matrixop()
{
	row = 0;
	col = 0;
	p = NULL;
}


matrixop :: ~matrixop()
{
	if ( p != NULL )
	delete [] p;
}

matrixop :: matrixop(const matrixop & mat)
{
int i;

row = mat.row;
col = mat.col;

p = new double[row*col];

for(i=0;i< ( mat.row )*( mat.col );i++)
p[i] = mat.p[i];

}


matrixop & matrixop :: operator=(const matrixop & mat)             //Assignment operator
{
int i;
row = mat.row;
col = mat.col;
k = mat.k;

if( p!= NULL )
delete [] p;

p = new double[row*col];

for(i=0;i<row*col;i++)
p[i] = mat.p[i];

return *this;
}


matrixop :: matrixop(int m, int n)
{

	p = new double[m*n];
	row = m;
	col = n;
	k = 0;
}

int matrixop :: r()
{
return row;
}

int matrixop :: c()
{
return col;
}


void matrixop :: disp()
{

int i,j;

for (i=0;i<row;i++)
{
for(j=0;j<col;j++)
printf("%10.8f ",p[j*row+i]);

cout << "\n";
}

}


void matrixop :: insert(double x)
{
p[k++] = x;
}

double matrixop :: get(int i)
{
return p[i];
}

void matrixop :: change(int i, double num)
{
p[i] = num;
}

void matrixop :: mul(matrixop & matA,char transA,matrixop & matB,char transB, double alpha)
{
int lda,ldb,ldc,beta;
int m,n,k;
if ( transA == 'N' && transB == 'N' )
{
(*this) = matrixop(matA.row,matB.col);
m = matA.row;
n = matB.col;
k = matA.col;
}
else if ( transA == 'N' && transB == 'T' )
{
(*this) = matrixop(matA.row,matB.row);
m = matA.row;
n = matB.row;
k = matA.col;
}
else if ( transA == 'T' && transB == 'N' )
{
(*this) = matrixop(matA.col,matB.col);
m = matA.col;
n = matB.col;
k = matA.row;
}
else if ( transA == 'T' && transB == 'T' )
{
(*this) = matrixop(matA.col,matB.row);
m = matA.col;
n = matB.row;
k = matA.row;
}

lda = matA.row;
ldb = matB.row;
ldc = row;
beta = 0;
dgemm_(&transA, &transB, &m, &n, &k, &alpha, matA.p, &lda, matB.p, &ldb, &beta, p, &ldc);

}

/*
matrixop & matrixop :: operator*(const matrixop & mat)
{
double alpha;
int lda,ldb,ldc,beta;
lda = row;
ldb = mat.row;
ldc = row;
alpha = 1;
beta = 0;
char transA = 'N';
char transB = 'N';
int M = row;
int N = mat.col;
int K = col;

matrixop *C = new matrixop;
(*C) = matrixop(row,mat.col);
dgemm_(&transA, &transB, &M, &N, &K, &alpha, p, &lda, mat.p, &ldb, &beta, (*C).p, &ldc);

return *C;
}
*/


double matrixop :: rcond()
{
int i,j;
double *array = new double[row*col];

if ( row != col )
cout << " NOT A SQUARE matrixop \n";

//array = new double[row*col];
for(i=0;i<row*col;i++)
array[i] = p[i];

int lda,INFO;
char NORM = '1';
double ANORM;
double RCOND,WORK[4*row];
int IWORK[row];

lda = row;
int IPIV[row];
double sum;

for (i=0;i<col;i++)
{
sum = 0;
for (j=0;j<row;j++)
sum = sum + p[col*i + j];

if ( ANORM < sum )
ANORM = sum;
}



dgetrf_(&row, &col, array, &lda, IPIV, &INFO);
dgecon_(&NORM, &row, array, &lda, &ANORM, &RCOND, WORK,IWORK,&INFO);



delete [] array;
return RCOND;
}

void matrixop :: inv()
{
int INFO,lda = row;
int IPIV[row+1];
int LWORK=row*row;
double WORK[row*row];

dgetrf_(&row, &col, p, &lda, IPIV, &INFO);

if (INFO == 0)
dgetri_(&row, p, &lda, IPIV, WORK, &LWORK, &INFO );
//else
//cout << "SOMETHING WRONG IN INVERSE";

}

void matrixop :: div(matrixop & mat)
{
int i;
double * array = new double[mat.row*mat.col];

for(i=0;i<mat.row*mat.col;i++)
array[i] = mat.p[i];


int NRHS = col;
int lda = mat.row;
int ldb = row;
int *IPIV = new int[row];
int INFO;

dgesv_(&mat.row,&NRHS,array,&lda,IPIV,p,&ldb,&INFO );

delete [] array;

delete [] IPIV;
}

double matrixop :: det()
{
int M,N,i;
int lda = row;
M = row;
N = col;
double * array = new double[row*col];
int IPIV[row];
int INFO;
double prod=1;

for(i=0;i<row*col;i++)
array[i] = p[i];

dgetrf_( &M, &N, array, &lda, IPIV, &INFO );

for(i=0;i<row*col;i=i+row+1)
prod = prod*array[i];

delete [] array;

return prod;
}


void matrixop :: pinv()
{
int i,j;
double * A = new double[row*col];
double tol = 1.8829e-013;
char JOBU = 'A', JOBVT = 'A';
int LDA,LDU,LWORK,INFO,LDVT;
int M = row;
int N = col;
LDA = row;
double * S = new double [ row <= col ? row : col ];
double * U, * VT, * WORK;
U = new double [row*col];
VT = new double [row*col];

for(i=0;i<row*col;i++)
A[i] = p[i];

LDA = row;
LDVT = col;
LDU = row;
LWORK = 5*row;
WORK = new double [LWORK];

dgesvd_( &JOBU, &JOBVT, &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &INFO );

matrixop Utran(row,col),VSinv(row,col);

for(i=0;i<col;i++)
{
	for(j=0;j<row;j++)
	{
		VSinv.insert(VT[j*col+i]*( S[i] > tol ? 1/S[i] : 0 ) );
		Utran.insert(U[j*col+i]);
	}
}

(*this).mul(VSinv,'N',Utran,'N',1);

delete [] S;
delete [] U;
delete [] VT;
delete [] WORK;
delete [] A;

}

#endif
