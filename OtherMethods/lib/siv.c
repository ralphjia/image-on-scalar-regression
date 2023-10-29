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
     #define RImg_IN prhs[0] /* Residual Image data n by N */
     #define XYZ_IN prhs[1] /* Coordinate of Image data N by 3 */
     #define H_IN prhs[2] /* bandwidth candidates NH by 1 */
     #define SRImg_OUT plhs[0] /* Smoothed Residual Image data n by N */
     #define h_OUT plhs[1] /* optimal bandwith scalar */ 
     #define GCV_OUT plhs[2] /* all GCVs */ 
 
	 double *matrixInRImg, *matrixInXYZ, *vectorInH, *matrixOutSRImg,*matrixOutGCV; 
            /* do not need to claim scalar output or creat pointer */
     double **RImg,**XYZ,**SRImg;
     double *H, *GCV;
     double minh;
     int n, N, NH;
     int i, j, k, l, ii;
     double **Zh, **KZh, **KZh2, **S, **ISR;
     double trS, canh, tkzh,minGCV, currentGCV;
     double **InvKZh2,*col;
     int *indx;    

	 matrixInRImg = mxGetPr(RImg_IN); /* n by N matrix */
     matrixInXYZ = mxGetPr(XYZ_IN); /* N by 3 matrix */   
     vectorInH = mxGetPr(H_IN); /* NH by 1 vector */ 
	
     n = mxGetM(RImg_IN); /* number of observations */
     N = mxGetN(RImg_IN); /* number of voxels */
     NH = mxGetNumberOfElements(H_IN); /* number of bandwidth */
	
	 SRImg_OUT = mxCreateDoubleMatrix(n,N,mxREAL);
     matrixOutSRImg = mxGetPr(SRImg_OUT); 
    
     GCV_OUT=mxCreateDoubleMatrix(1,NH,mxREAL);
     matrixOutGCV = mxGetPr(GCV_OUT);
          
          
     /* convert Matlab matrix to C matrix */
     RImg = dmatrix(1,n,1,N); /* Residaul image in C double pointer format */
     XYZ = dmatrix(1,N,1,3); /* Coordinate of voxels in C double pointer format */
     H = dvector(1,NH); /* Bandwidth Candidates in C single pointer format */
     SRImg = dmatrix(1,n,1,N); /* Smoothed Residaul image in C double pointer format */  
     GCV = dvector(1,NH); /* GCV in C double vector format */
     
     for (i=0;i<n;i++){
		for (j=0;j<N;j++) RImg[i+1][j+1]=matrixInRImg[i+j*n];
     }
     
     for (i=0;i<N;i++){
		for (j=0;j<3;j++) XYZ[i+1][j+1]=matrixInXYZ[i+j*N];
     } 
     
     for (i=0;i<NH;i++){
         H[i+1]=vectorInH[i];
     }
     /* convert Matlab matrix to C matrix */
    

     /* MAIN CODE --- Find minmum GCV and its h */
     
     Zh = dmatrix(1,4,1,N); 
     KZh = dmatrix(1,4,1,N); 
     KZh2 = dmatrix(1,4,1,4);
     S = dmatrix(1,N,1,N);
     ISR = dmatrix(1,N,1,n);
     
     
     InvKZh2 = dmatrix(1,4,1,4);
     col = dvector(1,4);
     indx = ivector(1,4);
     
     minGCV=10000000000000000.0;
     minh=0.0;
     for (i=1;i<=NH;i++){
         canh = H[i];
     
         /* calculate smooth matrix S */   
         for (j=1;j<=N;j++){
             /* calculate Zh and KZh matrices */
             for (k=1;k<=N;k++){
                 Zh[1][k]=1;
                 Zh[2][k]=(XYZ[k][1]-XYZ[j][1])/canh;
                 Zh[3][k]=(XYZ[k][2]-XYZ[j][2])/canh;
                 Zh[4][k]=(XYZ[k][3]-XYZ[j][3])/canh;
                 tkzh=exp(-0.5*(Zh[2][k]*Zh[2][k]+Zh[3][k]*Zh[3][k]+Zh[4][k]*Zh[4][k]))/canh;
                 KZh[1][k]=tkzh;
                 KZh[2][k]=tkzh*Zh[2][k];
                 KZh[3][k]=tkzh*Zh[3][k];
                 KZh[4][k]=tkzh*Zh[4][k];                 
             }         
             /* calculate Zh and KZh matrices */
         
             /* calculate KZhs2  matrix */       
             for (k=1;k<=4;k++){
                 for (l=1;l<=4;l++){
                     KZh2[k][l]=0.0;
                      for (ii=1;ii<=N;ii++){
                          KZh2[k][l]+=KZh[k][ii]*Zh[l][ii];
                      }    
                 }
             }
             /* calculate KZh2 matrix */
         
             /* find the inverse of KZh2, denoted by InvKZh2 */     
             ludcmp(KZh2,4,indx); 
             for (k=1;k<=4;k++) { 
                 for (l=1;l<=4;l++) {
                     col[l]=0.0;
                 }
                 col[k]=1.0;
                 lubksb(KZh2,4,indx,col);
                 for(l=1;l<=4;l++) {
                     InvKZh2[l][k]=col[l];
                 }
              }
            /* find the inverse of KZh2, denoted by InvKZh2 */ 
         
            for (k=1;k<=N;k++) {
                S[j][k]=0.0;
                for (l=1;l<=4;l++) { 
                    S[j][k]+=InvKZh2[1][l]*KZh[l][k];
                } 
            }
         }       
         /* calculate smooth matrix S */
         
          /* calculate trace trS */
         trS=0.0;
         for (j=1;j<=N;j++){
             trS+=S[j][j];
         }
         /* calculate trace trS */
         
         /* calculate matrix ISR */
         for (j=1;j<=N;j++){
             for (k=1;k<=n;k++) {
                 ISR[j][k]=RImg[k][j];
                 for (l=1;l<=N;l++) {
                     ISR[j][k]-=S[j][l]*RImg[k][l];
                 }    
             }
         }         
         /* calculate matrix ISR */
         
         /* calculate GCV */
         currentGCV=0.0;
         for (j=1;j<=n;j++){
             for (k=1;k<=N;k++){
                 currentGCV+=ISR[k][j]*ISR[k][j];
             }
         } 
         currentGCV=currentGCV/(1-trS/N)/(1-trS/N);
         /* calculate GCV */
         
         GCV[i]=currentGCV;
         if (currentGCV<minGCV) {
            minGCV=currentGCV;
            minh=canh;
         }
     }
     /* MAIN CODE  --- Find minmum GCV and its h */  
    
     /* MAIN CODE  --- calculate smoothed residual image */ 
         canh=minh;
         
         /* calculate smooth matrix S */   
         for (j=1;j<=N;j++){
             /* calculate Zh and KZh matrices */
             for (k=1;k<=N;k++){
                 Zh[1][k]=1;
                 Zh[2][k]=(XYZ[k][1]-XYZ[j][1])/canh;
                 Zh[3][k]=(XYZ[k][2]-XYZ[j][2])/canh;
                 Zh[4][k]=(XYZ[k][3]-XYZ[j][3])/canh;
                 tkzh=exp(-0.5*(Zh[2][k]*Zh[2][k]+Zh[3][k]*Zh[3][k]+Zh[4][k]*Zh[4][k]))/canh;
                 KZh[1][k]=tkzh;
                 KZh[2][k]=tkzh*Zh[2][k];
                 KZh[3][k]=tkzh*Zh[3][k];
                 KZh[4][k]=tkzh*Zh[4][k];                 
             }         
             /* calculate Zh and KZh matrices */
         
             /* calculate KZhs2  matrix */       
             for (k=1;k<=4;k++){
                 for (l=1;l<=4;l++){
                     KZh2[k][l]=0.0;
                      for (ii=1;ii<=N;ii++){
                          KZh2[k][l]+=KZh[k][ii]*Zh[l][ii];
                      }    
                 }
             }
             /* calculate KZh2 matrix */
         
             /* find the inverse of KZh2, denoted by InvKZh2 */     
             ludcmp(KZh2,4,indx); 
             for (k=1;k<=4;k++) { 
                 for (l=1;l<=4;l++) {
                     col[l]=0.0;
                 }
                 col[k]=1.0;
                 lubksb(KZh2,4,indx,col);
                 for(l=1;l<=4;l++) {
                     InvKZh2[l][k]=col[l];
                 }
              }
            /* find the inverse of KZh2, denoted by InvKZh2 */ 
         
            for (k=1;k<=N;k++) {
                S[j][k]=0.0;
                for (l=1;l<=4;l++) { 
                    S[j][k]+=InvKZh2[1][l]*KZh[l][k];
                } 
            }
         }       
         /* calculate smooth matrix S */
         
         
         /* calculate matrix SRImg */
         for (j=1;j<=n;j++){
             for (k=1;k<=N;k++) {
                 SRImg[j][k]=0.0;
                 for (l=1;l<=N;l++) {
                     SRImg[j][k]+=S[k][l]*RImg[j][l];
                 }    
             }
         }
         /* calculate matrix SRImg */     
     /* MAIN CODE  --- calculate smoothed residual image */
  
/*/    for (i=1;i<=n;i++){
//		for (j=1;j<=N;j++){ 
//            SRImg[i][j] = RImg[i][j];
//        }
//    }*/
   
         
    /* transfer output to matlab form */
    for (i=0;i<n;i++){
		for (j=0;j<N;j++){ 
            matrixOutSRImg[i+j*n] = SRImg[i+1][j+1]; 
        }
	}
    
     h_OUT = mxCreateDoubleScalar(minh);
     

    for (i=0;i<NH;i++){
	    matrixOutGCV[i] = GCV[i+1]; 
	}
     
    /* transfer output to matlab form */
    
    free_dmatrix(RImg,1,n,1,N); 
    free_dmatrix(XYZ,1,N,1,3);
    free_dvector(H,1,NH);
    free_dvector(GCV,1,NH);
    free_dmatrix(SRImg,1,n,1,N); 
    free_dmatrix(Zh,1,4,1,N); 
    free_dmatrix(KZh,1,4,1,N);
    free_dmatrix(KZh2,1,4,1,4);
    free_dmatrix(S,1,N,1,N);
    free_dmatrix(ISR,1,N,1,n);
    free_dmatrix(InvKZh2,1,4,1,4);
    free_dvector(col,1,4);
    free_ivector(indx,1,4);    
    
	return;
}
