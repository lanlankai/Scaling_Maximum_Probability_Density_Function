/*
* Yongxiang Huang, last modification: 23/03/2010
* yongxianghuang@gmail.com
*
*/

/* This function is to estimate the pdf scaling of velocity increments*/

#include <stdlib.h>
#include <stdio.h>
#include "mex.h"
#include <math.h>



/************************************************************************/
/*                                                                      */
/* MAIN FUNCTION                                                        */
/*                                                                      */
/************************************************************************/

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
  
    /* declarations */
    int i,j,ntmp,n,Nx,Ntau,nbin;
    
    double *x,*pdf,*tau,tmp1,tmp2,*tmp3,*tmpdf,*bin; /*given time seris, time delay, maximum q, dq is the step*/
      
    
/*     check input*/
    if (nrhs!=2)     mexErrMsgTxt("You have to input four parameters!");
    if (mxIsEmpty(prhs[0]))mexErrMsgTxt("Time SERIES is empty!");
    if (mxIsEmpty(prhs[1]))mexErrMsgTxt("Time Delay is empty!");
  
    /* get input data */
    x=mxGetPr(prhs[0]);
    tau=mxGetPr(prhs[1]);
    
    Nx=mxGetN(prhs[0]);
    ntmp=mxGetM(prhs[0]);
    if (ntmp>Nx) Nx=ntmp;
    
    
    Ntau=mxGetN(prhs[1]);
    ntmp=mxGetM(prhs[1]);
    if (ntmp>Ntau) Ntau=ntmp;
    
    /*specify bins*/
    nbin=22;
    bin=(double *)malloc(nbin*sizeof(double)); /*specify  bins*/
    for(i=0;i<nbin;i++)
    {
        bin[i]=-0.525+0.05*i;
    }
    nbin=nbin-1;
    
    tmpdf=(double *)malloc(nbin*sizeof(double)); /*specify  pdf*/
    
    tmp3=(double *)malloc(Nx*sizeof(double)); 
    
    
    plhs[0]=mxCreateDoubleMatrix(2,Ntau,mxREAL);
    pdf=mxGetPr(plhs[0]);
    for(i=0;i<2*Ntau;i++)pdf[i]=0.0;
    for(i=0;i<Ntau;i++)/*the bigest loop for each tau*/
  {
        
        for(j=0;j<nbin;j++)tmpdf[j]=0.0;
        ntmp=Nx-tau[i];
        n=tau[i];
        
        for(j=0;j<ntmp;j++)tmp3[j]=x[n+j]-x[j]; /*get velocity increments*/
        tmp2=0.0;
        for(j=0;j<ntmp;j++){tmp2=tmp2+tmp3[j]*tmp3[j];}
        tmp2=tmp2/ntmp;
        tmp2=pow(tmp2, 0.5);

        for(j=0;j<ntmp;j++)tmp3[j]=tmp3[j]/tmp2;
        
        for(j=0;j<ntmp;j++) /*loop for pdf estimation*/
        {
            for(n=0;n<nbin-1;n++)
            {
                if(tmp3[j]>=bin[n] && tmp3[j]<bin[n+1])
                {
                    tmpdf[n]=tmpdf[n]+1;
                    break;
                }
            
            }
                       
        }
        /*pdf[2*i]=tmpdf[0];*/
        for(j=0;j<nbin-1;j++) 
        {
            if (tmpdf[j]>pdf[2*i]) pdf[2*i]=tmpdf[j];
        }
        
        pdf[2*i]=pdf[2*i]/ntmp/tmp2;
        
        for(j=0;j<ntmp-1;j++)
        {
            if(tmp3[j]*tmp3[j+1]<0) {pdf[i*2+1]=pdf[i*2+1]+1.0;}
        }
        pdf[2*i+1]=pdf[2*i+1]/ntmp;
  }
    
    free(bin);
    free(tmp3);
    free(tmpdf);
}