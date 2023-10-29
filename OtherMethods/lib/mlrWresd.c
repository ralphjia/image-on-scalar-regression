#include <math.h>
#include "mex.h"
#define TINY 1.0e-30 /* A small number. */


double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	double **m;

	m=(double **) mxMalloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) mxMalloc((unsigned) (nch-ncl+1)*sizeof(double));
		m[i] -= ncl;
	}
	return m;
}

void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) mxFree((char*) (m[i]+ncl));
	mxFree((char*) (m+nrl));
}

double *dvector(int nl, int nh)
{
	double *v;

	v=(double *) mxMalloc((unsigned) (nh-nl+1)*sizeof(double));
	return v-nl;
}

void free_dvector(double *v, int nl, int nh)
{
	mxFree((char*) (v+nl));
}

int *ivector(int nl, int nh)
{
	int *v;

	v=(int *) mxMalloc((unsigned) (nh-nl+1)*sizeof(int));
	return v-nl;
}

void free_ivector(int *v, int nl, int nh)
{
	mxFree((char*) (v+nl));
}

void ludcmp(double **a, int n, int *indx)
{
int i,imax,j,k;
double big,dum,sum,temp;
double *vv; 

vv=dvector(1,n);
/*d=1.0; */
for (i=1;i<=n;i++) {
	big=0.0; 
    for (j=1;j<=n;j++)
        if ((temp=fabs(a[i][j])) > big) big=temp;
        vv[i]=1.0/big; 
}
for (j=1;j<=n;j++) { 
    for (i=1;i<j;i++) { 
    sum=a[i][j];
    for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
    a[i][j]=sum;
}
big=0.0; 
for (i=j;i<=n;i++) { 
     sum=a[i][j]; 
     for (k=1;k<j;k++)
         sum -= a[i][k]*a[k][j];
     a[i][j]=sum;
     if ( (dum=vv[i]*fabs(sum)) >= big) {
        big=dum;
        imax=i;
     }
}
if (j != imax) { 
   for (k=1;k<=n;k++) { 
       dum=a[imax][k];
       a[imax][k]=a[j][k];
       a[j][k]=dum;
   }
  /*d = -(*d); */
  vv[imax]=vv[j]; 
}
indx[j]=imax;
if (a[j][j] == 0.0) a[j][j]=TINY;
if (j != n) { 
   dum=1.0/(a[j][j]);
   for (i=j+1;i<=n;i++) a[i][j] *= dum;
}
} 
free_dvector(vv,1,n);
}

void lubksb(double **a, int n, int *indx, double *b)
{
int i,ii=0,ip,j;
double sum;
for (i=1;i<=n;i++) { 
ip=indx[i];
sum=b[ip];
b[ip]=b[i];
if (ii>0)
    for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum>0) ii=i; 
	b[i]=sum; 
}
for (i=n;i>=1;i--) { 
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i]; 
} 
}


/*void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]){*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
     #define X_IN prhs[0] /* design matrix */
     #define Y_IN prhs[1] /* response matrix */
     #define B_OUT plhs[0] /* parameter matrix */
     #define C_OUT plhs[1] /* covariances of parameters */
     #define R_OUT plhs[2] /* residuals of response */
     
	 double *matrixInX, *matrixOutB, *matrixInY, *matrixOutC, *matrixOutR;
     double **XX,**X,**Y; 
     int n, p, N;
	 int i, j, k, l;
     double **IXX,*col;
     int *indx; 

	 matrixInX = mxGetPr(X_IN); /* n by p matrix */
     matrixInY = mxGetPr(Y_IN); /* n by N matrix */   
	
     n = mxGetM(X_IN); /* number of observations */
     p = mxGetN(X_IN); /* number of covariates */
     N = mxGetN(Y_IN); /* number of voxels */
	
	 B_OUT = mxCreateDoubleMatrix(p,N,mxREAL);
     C_OUT = mxCreateDoubleMatrix(p,N,mxREAL);
     R_OUT = mxCreateDoubleMatrix(n,N,mxREAL);
     
     matrixOutB = mxGetPr(B_OUT); 
     matrixOutC = mxGetPr(C_OUT);
     matrixOutR = mxGetPr(R_OUT);
    
     /* convert Matlab matrix to C matrix */
     X = dmatrix(1,n,1,p); /* design matrix in C double pointer format */
     Y = dmatrix(1,n,1,N); /* response matrix in C double pointer format */
     
     for (i=0;i<n;i++){
		for (j=0;j<p;j++) X[i+1][j+1]=matrixInX[i+j*n];
     }
     
     for (i=0;i<n;i++){
		for (j=0;j<N;j++) Y[i+1][j+1]=matrixInY[i+j*n];
     }     
     /* convert Matlab matrix to C matrix */
    
     /* calculate XX matrix */
     XX = dmatrix(1,p,1,p);           
     for (i=1;i<=p;i++){
		for (j=1;j<=p;j++){
            XX[i][j] = 0.0; 
            for (k=1;k<=n;k++){
                 XX[i][j] += X[k][i]*X[k][j];
            }
        }
	 }      
    /* calculate XX matrix */
     
     /* find the inverse of XX, denoted by IXX */     
    
     IXX = dmatrix(1,p,1,p);
     col = dvector(1,p);
     indx = ivector(1,p);
          
     ludcmp(XX,p,indx); 
     for (j=1;j<=p;j++) { 
         for (i=1;i<=p;i++) col[i]=0.0;
              col[j]=1.0;
              lubksb(XX,p,indx,col);
              for(i=1;i<=p;i++) IXX[i][j]=col[i];
      }
      /* find the inverse of XX, denoted by IXX */  
     
        
    /* calculate projection matrix PX and hessian matrix HX */     
    double **PX, **HX; /* note here, HX = I-H = I-X(X'X)^{-1}X' */
    PX = dmatrix(1,p,1,n);
    HX = dmatrix(1,n,1,n);
    
    for (i=1;i<=p;i++){
        for (j=1;j<=n;j++){
            PX[i][j] = 0.0;
            for (k=1;k<=p;k++) PX[i][j] +=IXX[i][k]*X[j][k]; 
        }
    }
     
    for (i=1;i<=n;i++){
        for (j=1;j<=n;j++){
            HX[i][j] = 0.0;
            for (k=1;k<=p;k++) HX[i][j] -=X[i][k]*PX[k][j]; /* -X(X'X)^{-1}X' */
            if (i == j) HX[i][j]++; /* I-X(X'X)^{-1}X' */
        }
    }
    /* calculate projection matrix PX and hessian matrix HX */ 
    
    /* calculate parameter Beta, its covariance Covb and Residuals R */
    double **Beta, **Covb, **R;
    double sigma2,tempR;
    Beta = dmatrix(1,p,1,N);
    Covb = dmatrix(1,p,1,N);
    R = dmatrix(1,n,1,N);    
    
    for (i=1;i<=p;i++){
        for (j=1;j<=N;j++){
            Beta[i][j] = 0.0;
            for (k=1;k<=n;k++) Beta[i][j] += PX[i][k]*Y[k][j];
        }
    }
    
    
    for (i=1;i<=N;i++){
        sigma2 = 0.0;
        for (j=1;j<=n;j++){
             for (k=1;k<=n;k++){   
		         sigma2 += Y[j][i]*HX[j][k]*Y[k][i];
             }
        }
        sigma2 = sigma2/(n-p);
        for (l=1;l<=p;l++) Covb[l][i] = sigma2*IXX[l][l];           
    } 
    
    for (i=1;i<=N;i++){
        for (j=1;j<=n;j++){
            tempR = 0.0;
            for (k=1;k<=n;k++){   
		        tempR += HX[j][k]*Y[k][i];
            }
            R[j][i] = tempR;
       }
    }
     /* calculate parameter Beta, its covariance Covb and Residuals R */
    
    /* transfer output to matlab form */
    for (i=0;i<p;i++){
		for (j=0;j<N;j++){ 
            matrixOutB[i+j*p] = Beta[i+1][j+1]; 
        }
	}
    
    for (i=0;i<p;i++){
		for (j=0;j<N;j++){ 
            matrixOutC[i+j*p] = Covb[i+1][j+1]; 
        }
	}    
    
    for (i=0;i<n;i++){
		for (j=0;j<N;j++){ 
            matrixOutR[i+j*n] = R[i+1][j+1]; 
        }
	}
    /* transfer output to matlab form */
    
    free_dmatrix(X,1,n,1,p); 
    free_dmatrix(Y,1,n,1,N);
   	free_dmatrix(XX,1,p,1,p);
    free_dmatrix(IXX,1,p,1,p);
    free_dvector(col,1,p);
    free_ivector(indx,1,p);
    free_dmatrix(PX,1,p,1,n);
    free_dmatrix(HX,1,n,1,n);
    free_dmatrix(Beta,1,p,1,N);
    free_dmatrix(Covb,1,p,1,N);
    free_dmatrix(R,1,n,1,N);
    
	return;
}




