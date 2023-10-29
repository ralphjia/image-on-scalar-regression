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

int **imatrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	int **m;

	m=(int **) mxMalloc((unsigned) (nrh-nrl+1)*sizeof(int*));
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(int *) mxMalloc((unsigned) (nch-ncl+1)*sizeof(int));
		m[i] -= ncl;
	}
	return m;
}

void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) mxFree((char*) (m[i]+ncl));
	mxFree((char*) (m+nrl));
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

/********************************************************************************************/
/*/ revised @ June 16, 2011 */
/*/ correct some errors on output coverting */
/*/ compact transfering between C and Matlab format*/
/*/ redecued two output entries*/
/*/ revised @ June 17, 2011*/
/*/ add one more input on chi-squared critical values*/
/*/ revised @ July 21, 2011*/
/*/ add error variance in the adaptive smoothing'*/
/*/ fixed error that did not add sigma2_epsilon */
/* at the first step of calculating sigma2Last @ March 12, 2012 */
/*********************************************************************************************/


/*void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]){*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
     #define BETA_IN prhs[0] /* Initial Beta, N by p */
     #define ETA_IN prhs[1] /* Smoothed initial residaul, N by n */
     #define HS_IN prhs[2] /*  Adaptive radius candidate, nh by 1 */
     #define VNS_IN prhs[3] /* Number of maximal nearby voxel, N by 1  */
     #define VNSID_IN prhs[4] /*  Nearby voxel sequence ID, N by m */
     #define VDS_IN prhs[5] /*  Nearby voxel distance sequence, N by m */
     #define X_IN prhs[6] /*  covariate matrix, n by p */  
     #define CHISEQ_IN prhs[7] /*  chi critial value sequence, nh-2 by 1 */    
     #define SE_IN prhs[8] /*  error covariance, N by 1 */  
     #define CN_IN prhs[9] /*  a scalar */  
     #define BETA_OUT plhs[0] /* Final Smoothed Beta, N by p */
     #define HS_OUT plhs[1] /* Radius for each voxel, N by p */ 
     #define TEST_OUT plhs[2] 
 
	 double *matrixInBeta, *matrixInEta, *vectorInHs, *vectorInVns, *matrixInVnsid, 
            *matrixInVds, *matrixInX, *vectorInChiseq, *vectorInSe, *scalarInCn, 
            *matrixOutBeta, *matrixOutHs, *matrixOutTest; 
     double **BetaIn, **EtaIn, **VnsidIn, **VdsIn, **XIn, *SeIn, CnIn,
            **BetaOut, **HsOut;
     double *HsIn, *VnsIn, *ChiseqIn;
     double **XX; 
     double **IXX,*col;
     int *indx;  
     int n, N, p, nh, m, nh1,dh;
     int i, j, k, l, ii, jj;
     double **BetaFix, **Sigma2Fix, **BetaLast, **Sigma2Last, **Sigma2Out, **wBeta;
     double *sumwBeta;
     double tSigma2, dBeta, Kloc, Kst, twBeta, Cn, tDist;
     int **stopUpdate;
     int nbycount, tId1, tId2;

	 matrixInBeta = mxGetPr(BETA_IN); /* N by p matrix */
     matrixInEta = mxGetPr(ETA_IN); /* N by n matrix */   
     vectorInHs = mxGetPr(HS_IN); /* nh by 1 vector */ 
     vectorInVns = mxGetPr(VNS_IN); /* N by 1 matrix */
     matrixInVnsid = mxGetPr(VNSID_IN); /* N by m matrix */
     matrixInVds = mxGetPr(VDS_IN); /* N by m matrix */
     matrixInX = mxGetPr(X_IN); /* n by p matrix */
     vectorInChiseq = mxGetPr(CHISEQ_IN); /* nh-dh by 1 vector */
     vectorInSe = mxGetPr(SE_IN); /* N by 1 vector */
     scalarInCn = mxGetPr(CN_IN); /**/
	
     n = mxGetN(ETA_IN); /* number of observations */
     N = mxGetM(ETA_IN); /* number of voxels */
     nh = mxGetNumberOfElements(HS_IN); /* number of radiu candidate */
     nh1 = mxGetNumberOfElements(CHISEQ_IN); /* number of check points */
     p = mxGetN(BETA_IN); /* number of covariates */
     m = mxGetN(VNSID_IN); /* number of nearby voxels */
	
	 BETA_OUT = mxCreateDoubleMatrix(N,p,mxREAL);
     HS_OUT = mxCreateDoubleMatrix(N,p,mxREAL); 
     
     matrixOutBeta = mxGetPr(BETA_OUT);
     matrixOutHs = mxGetPr(HS_OUT);  
     
     TEST_OUT = mxCreateDoubleMatrix(N,p,mxREAL);
     matrixOutTest = mxGetPr(TEST_OUT);   
     
          
     /* convert Matlab matrix to C matrix */
     BetaIn = dmatrix(1,N,1,p); /* Parameter initial map in C double pointer format */
     EtaIn = dmatrix(1,N,1,n); /* Smoothed initial residaul in C double pointer format */
     HsIn = dvector(1,nh); /* Adaptive radius candidate in C single pointer format */
     VnsIn = dvector(1,N); /* Number of maximal nearby voxel in C single pointer format */
     VnsidIn = dmatrix(1,N,1,m); /* Nearby voxel sequence ID in C double pointer format */  
     VdsIn = dmatrix(1,N,1,m); /* Nearby voxel distance sequence in C double pointer format */  
     XIn = dmatrix(1,n,1,p); /* Covariate matrix in C double pointer format */ 
     ChiseqIn = dvector(1,nh1); /* Chi square critical vector in C double pointer format */ 
     SeIn = dvector(1,N); /* Error covariance vector in C double pointer format */ 
     
     for (i=0;i<N;i++){
		for (j=0;j<p;j++) BetaIn[i+1][j+1]=matrixInBeta[i+j*N];
     }
     
     for (i=0;i<N;i++){
		for (j=0;j<n;j++) EtaIn[i+1][j+1]=matrixInEta[i+j*N];
     } 
     
     for (i=0;i<nh;i++){
         HsIn[i+1]=vectorInHs[i];
     }
     
     for (i=0;i<N;i++){
         VnsIn[i+1]=vectorInVns[i];
     }
     
     for (i=0;i<N;i++){
		for (j=0;j<m;j++){ 
            VnsidIn[i+1][j+1]=matrixInVnsid[i+j*N];
            VdsIn[i+1][j+1]=matrixInVds[i+j*N];
        }
     }
     
     for (i=0;i<n;i++){
		for (j=0;j<p;j++) XIn[i+1][j+1]=matrixInX[i+j*n];
     }
     
     for (i=0;i<nh1;i++){
         ChiseqIn[i+1]=vectorInChiseq[i];
     }
     
     for (i=0;i<N;i++){
         SeIn[i+1]=vectorInSe[i];
     }
     
     CnIn = scalarInCn[0];
     
     /* convert Matlab matrix to C matrix */
    

     /* --- MAIN CODE ---  */
     BetaOut = dmatrix(1,N,1,p); /* Parameter final map in C double pointer format */
     HsOut = dmatrix(1,N,1,p); /* Final radius map in C double pointer format */
             
     /* calculate XX matrix */
       
     XX = dmatrix(1,p,1,p);           
     for (i=1;i<=p;i++){
		for (j=1;j<=p;j++){
            XX[i][j] = 0.0; 
            for (k=1;k<=n;k++){
                 XX[i][j] += XIn[k][i]*XIn[k][j];
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
     
     
     BetaFix = dmatrix(1,N,1,p);
     Sigma2Fix = dmatrix(1,N,1,p);
     BetaLast = dmatrix(1,N,1,p);
     Sigma2Last = dmatrix(1,N,1,p);
     Sigma2Out = dmatrix(1,N,1,p);
     wBeta = dmatrix(1,m,1,p);
     sumwBeta = dvector(1,p);
     stopUpdate = imatrix(1,N,1,p);
     dh = nh-nh1;

         
     /* find BetaFix and Sigma2Fix at step 3 */
     for (i=1;i<=dh;i++){

         /* initialize BetaLast and Sigma2Last */
         if (i==1)
         {
             for (j=1;j<=N;j++)
             {
                 tSigma2 = 0.0;
                 for (k=1;k<=n;k++)
                 {
                      tSigma2 += EtaIn[j][k]*EtaIn[j][k];
                 }
                 tSigma2 = tSigma2/(n-p) + SeIn[j];
                 for (k=1;k<=p;k++)
                 {
                     BetaLast[j][k] = BetaIn[j][k];
                     Sigma2Last[j][k] = tSigma2*IXX[k][k];                     
                 }
             }
         } 
         else 
         {
              for (j=1;j<=N;j++)
              {
                 for (k=1;k<=p;k++)
                 {
                     BetaLast[j][k] = BetaOut[j][k];
                     Sigma2Last[j][k] = Sigma2Out[j][k];                     
                 }
             }
         }
         /* initialize BetaLast and Sigma2Last */
         
         /* iterate all voxels to obtain BetaOut and Sigma2Out */
         for (j=1;j<=N;j++){
         
             /* find nearby voxel*/
             nbycount = 0;
             for (k=1;k<=VnsIn[j];k++)
             {                 
                 if (VdsIn[j][nbycount+1] <= HsIn[i])
                 {
                     nbycount++;
                 }
             }
             /* find nearby voxel*/
             
             /* find weights of nearby voxel*/
             for (k=1;k<=p;k++)
             {
                 sumwBeta[k] = 0.0;
             }
                         
             for (k=1;k<=nbycount;k++)
             {
                 Kloc = 1-VdsIn[j][k]/HsIn[i]; /* negative?? */
                 tId1 = VnsidIn[j][k];
                 for (l=1;l<=p;l++)
                 {
                     dBeta = (BetaLast[j][l]-BetaLast[tId1][l])*(BetaLast[j][l]-BetaLast[tId1][l])/Sigma2Last[j][l]/(CnIn/(11-i));  
                     Kst = exp(-dBeta);
/*/                      Kst = 2.0*(1.0-dBeta*dBeta);*/
/*/                      if (Kst <= 0.0)*/
/*/                      {*/
/*/                          Kst = 0.0;*/
/*/                      }else if (Kst >= 1.0)*/
/*/                      {*/
/*/                          Kst = 1.0;*/
/*/                      }else*/
/*/                      {*/
/*/                          Kst = Kst;*/
/*/                      }*/
                     wBeta[k][l] = Kloc*Kst;  
                     sumwBeta[l] = sumwBeta[l]+wBeta[k][l];  
                 }
             }             
             /* find weights of nearby voxel*/
             
             /* update BetaOut and Sigma2Out */
             for (k=1;k<=p;k++)
             {
                 BetaOut[j][k] = 0.0;
                 Sigma2Out[j][k] = 0.0;
                 for (l=1;l<=nbycount;l++)
                 {
                     tId1 = VnsidIn[j][l];
                     BetaOut[j][k] = BetaOut[j][k] + wBeta[l][k]*BetaLast[tId1][k];
                     for (ii=1;ii<=nbycount;ii++)
                     {
                         tId2 = VnsidIn[j][ii];
                         tSigma2 = 0.0;
                         for (jj=1;jj<=n;jj++)
                         {
                             tSigma2 = tSigma2 + EtaIn[tId1][jj]*EtaIn[tId2][jj];
                         }
                         tSigma2 = tSigma2/(n-p)*IXX[k][k];
                         if (l==ii)
                         {
                             tSigma2 = tSigma2 + SeIn[l]*IXX[k][k];
                         }    
                         Sigma2Out[j][k] = Sigma2Out[j][k] + wBeta[l][k]*wBeta[ii][k]*tSigma2;
                     }
                 }
                 BetaOut[j][k] = BetaOut[j][k]/sumwBeta[k];
                 Sigma2Out[j][k] = Sigma2Out[j][k]/(sumwBeta[k]*sumwBeta[k]);
             } 
             /* update BetaOut and Sigma2Out */
             
         }
         /* iterate all voxels to obtain BetaOut and Sigma2Out*/
        
     }  
     /* find BetaFix and Sigma2Fix at step 3 */
     
     
     /* set up BetaFix and Sigma2Fix, and stopUpdate */
     for (i=1;i<=N;i++)
     {
         for (j=1;j<=p;j++)
         {
             BetaFix[i][j] = BetaOut[i][j];
             Sigma2Fix[i][j] = Sigma2Out[i][j];
             stopUpdate[i][j] = 0;
         }
     }   
     /* set up BetaFix and Sigma2Fix, and stopUpdate */
     
     
     
     /* iterate to find the final Beta and Sigma2, and optimal h */
     for (i=dh+1;i<=nh;i++)
     {
         /* update BetaLast and Sigma2Last */
         for (j=1;j<=N;j++)
         {
             for (k=1;k<=p;k++)
             {
                 if (stopUpdate[j][k]==0) 
                 {    
                   BetaLast[j][k] = BetaOut[j][k];
                   Sigma2Last[j][k] = Sigma2Out[j][k];     
                 }
             }
         }
         /* update BetaLast and Sigma2Last */
         
         /* iterate all voxels to update BetaOut and Sigma2Out according to HsOut and stopUpdate */
         for (j=1;j<=N;j++){  
             
             /* find nearby voxel */
             nbycount = 0;
             for (k=1;k<=VnsIn[j];k++)
             {                 
                 if (VdsIn[j][nbycount+1] <= HsIn[i])
                 {
                     nbycount++;
                 }
             }
             /* find nearby voxel */
             
             /* find weights of nearby voxel */
             for (k=1;k<=p;k++)
             {
                 sumwBeta[k] = 0.0;
             }
                         
             for (k=1;k<=nbycount;k++)
             {
                 Kloc = 1-VdsIn[j][k]/HsIn[i];
                 tId1 = VnsidIn[j][k];
                 for (l=1;l<=p;l++)
                 {
                     dBeta = (BetaLast[j][l]-BetaLast[tId1][l])*(BetaLast[j][l]-BetaLast[tId1][l])/Sigma2Last[j][l]/(CnIn/(11-i));  
                     Kst = exp(-dBeta);                      
/*/                      Kst = 2.0*(1.0-dBeta*dBeta);*/
/*/                      if (Kst <= 0.0)*/
/*/                      {*/
/*/                          Kst = 0.0;*/
/*/                      }else if (Kst >= 1.0)*/
/*/                      {*/
/*/                          Kst = 1.0;*/
/*/                      }else*/
/*/                      {*/
/*/                          Kst = Kst;*/
/*/                      }*/
                     wBeta[k][l] = Kloc*Kst;  
                     sumwBeta[l] = sumwBeta[l]+wBeta[k][l];  
                 }
             }             
             /* find weights of nearby voxel */
             
             /* update BetaOut and Sigma2Out */
             for (k=1;k<=p;k++)
             {
                 if (stopUpdate[j][k]==0)
                 {
                     BetaOut[j][k] = 0.0;
                     Sigma2Out[j][k] = 0.0;
                     for (l=1;l<=nbycount;l++)
                     {
                         tId1 = VnsidIn[j][l];
                         BetaOut[j][k] = BetaOut[j][k] + wBeta[l][k]*BetaLast[tId1][k];
                         for (ii=1;ii<=nbycount;ii++)
                         {
                             tId2 = VnsidIn[j][ii];
                             tSigma2 = 0.0;
                             for (jj=1;jj<=n;jj++)
                             {
                                 tSigma2 = tSigma2 + EtaIn[tId1][jj]*EtaIn[tId2][jj];
                             }
                             tSigma2 = tSigma2/(n-p)*IXX[k][k];
                             if (l==ii)
                             {
                                tSigma2 = tSigma2 + SeIn[l]*IXX[k][k];
                             } 
                             Sigma2Out[j][k] = Sigma2Out[j][k] + wBeta[l][k]*wBeta[ii][k]*tSigma2;
                         }
                     }
                     BetaOut[j][k] = BetaOut[j][k]/sumwBeta[k];
                     Sigma2Out[j][k] = Sigma2Out[j][k]/(sumwBeta[k]*sumwBeta[k]);
                 }
/*/                  else if (stopUpdate[j][k]==1)*/
/*/                  {*/
/*/                      BetaOut[j][k] = BetaLast[j][k];*/
/*/                      Sigma2Out[j][k] = Sigma2Last[j][k];*/
/*/                  }                     */
             } 
             /* update BetaOut and Sigma2Out */
             
         }
        /* iterate all voxels to update BetaOut and Sigma2Out according to HsOut */
         
        /* update stopUpdate and HsOut */
         for (j=1;j<=N;j++)
         {
            for (k=1;k<=p;k++)
            {
                if (stopUpdate[j][k]==0)
                {
                     tDist = (BetaFix[j][k]-BetaOut[j][k])*(BetaFix[j][k]-BetaOut[j][k])/Sigma2Fix[j][k];
                    if (tDist > ChiseqIn[i-3])
                    {
                         stopUpdate[j][k] = 1;
                         HsOut[j][k] = HsIn[i];
                    }
                }                     
            }
         }
         /* update stopUpdate and HsOut */        
         
     }
     /* iterate to find the final Beta and Sigma2, and optimal h  and stopUpdate */
     
     
     /* update final stopUpdate and HsOut */
     for (j=1;j<=N;j++)
     {
         for (k=1;k<=p;k++)
         {
             if (stopUpdate[j][k]==0)
             {
                 stopUpdate[j][k] = 1;
                 HsOut[j][k] = HsIn[nh];
             }
         }                     
     }
    /* update final stopUpdate and HsOut */  
     
         /* update final BetaOut and Sigma2Out */
         for (j=1;j<=N;j++)
         {
             for (k=1;k<=p;k++)
             {
                 if (stopUpdate[j][k]==1) 
                 {    
                   BetaOut[j][k] = BetaLast[j][k];
                   Sigma2Out[j][k] = Sigma2Last[j][k];     
                 }
             }
         }
         /* update final BetaOut and Sigma2Out */
     
           
     /* --- MAIN CODE  --- */
     
      
    /* transfer output to matlab form */
    for (i=0;i<N;i++){
		for (j=0;j<p;j++){ 
            matrixOutBeta[i+j*N] = BetaOut[i+1][j+1]; 
            matrixOutHs[i+j*N] = HsOut[i+1][j+1]; 
        }
	}    
   
    for (i=0;i<N;i++){
		for (j=0;j<p;j++){ 
            matrixOutTest[i+j*N] = Sigma2Out[i+1][j+1]; 
        }
	}    
   
     /* transfer output to matlab form */
    
    
    free_dmatrix(BetaIn,1,N,1,p); 
    free_dmatrix(EtaIn,1,N,1,n);
    free_dvector(HsIn,1,nh);
    free_dvector(VnsIn,1,N);
    free_dmatrix(VnsidIn,1,N,1,m); 
    free_dmatrix(VdsIn,1,N,1,m); 
    free_dmatrix(XIn,1,n,1,p);
    free_dmatrix(BetaOut,1,N,1,p); 
    free_dmatrix(HsOut,1,N,1,p);  
    free_dvector(ChiseqIn,1,nh-3);
    free_dvector(SeIn,1,N);
    free_dmatrix(XX,1,p,1,p);
    free_dmatrix(IXX,1,p,1,p);
    free_dvector(col,1,p);
    free_ivector(indx,1,p);
    free_dmatrix(BetaFix,1,N,1,p);
    free_dmatrix(Sigma2Fix,1,N,1,p);
    free_dmatrix(BetaLast,1,N,1,p);
    free_dmatrix(Sigma2Last,1,N,1,p);
    free_dmatrix(Sigma2Out,1,N,1,p);
    free_dmatrix(wBeta,1,m,1,p);
    free_dvector(sumwBeta,1,p);   
    free_imatrix(stopUpdate,1,N,1,p);
   
	return;
}
