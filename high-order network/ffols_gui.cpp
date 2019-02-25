#include "mex.h"
#include "matrix.h"
double cal_dot(double *x,double* y, size_t N)
{
    double sums =0,csum1 =0,csum2 =0;
    int i =0;
    for (i =0;i<N;i++){
    sums  += x[i]*y[i];
    csum1 += x[i]*x[i];
    csum2 += y[i]*y[i];
    }
    if(csum1<0.01 || csum2<0.01)
        return 0;
    else
    return (sums*sums)/(csum1*csum2);
};
double ddot(double *x,double* y, size_t N)
{
    double sums =0;
    int i =0;
    for (i =0;i<N;i++){
    sums += x[i]*y[i];
    }
    return sums;
};
void ols(double *x,double *y, size_t N)
{
    double pw,wk;
    int i=0;
    pw = ddot(x,y,N);
    wk = ddot(x,x,N);
    if (wk>0.01)
    {       
        for(i=0;i<N;i++){
            y[i] = y[i]-(pw/wk)*x[i];
        }
    }
    return;
};
int minind(double *x,size_t N)
{        
    double i=x[0];
    int t=0;
    for(int j=0;j<N;j++)
    {
        if(i>x[j])
        {
            i = x[j];
            t = j;
        }
    }
    return t;
};

// double err_cal(double *Y, double *P,size_t M, size_t N, int kind, int * index,int *len,double *errs,double *pser)
double err_cal(double *Y, double *P,size_t M, size_t N, int kind, int * index,int *len,double *errs,double* threshold)
{

  //  double lambda = 3;
    double err = 0, maxerr=0,terr=0;
    int i =0 ,j = 0;
    int tk = 0;
    
//     int * nind = (int*)mxCalloc(N,sizeof(int));
    
    double *flags = (double*)mxCalloc(N,sizeof(double));
 
//     for (j=0;j<N;j++)  nind[j] = j;         
    for (i=0;i<N;i++)
    {
        err = cal_dot(Y+kind*M,P+i*M,M);
        if(maxerr<err)
        {
            maxerr = err;
            index[0] = i;
        }
        flags[i] = 0;
    }
    flags[index[0]] = 1;
    //pser[0] = (1-maxerr)/((1-lambda/M)*(1-lambda/M));
    err = maxerr;
    errs[0] = maxerr;
    for (i=1;i<N;i++)
    {
        maxerr = -1;
        for(j=0;j<N;j++)
        {
            if(flags[j]<0.01)
            {
                ols(P+(index[i-1])*M,P+j*M,M);
                terr = cal_dot(Y+kind*M,P+j*M,M);
                if (maxerr<terr)
                {
                    index[i] = j;
                    maxerr = terr;
					//mexPrintf("%d ",index[i]);
					//mexPrintf("%f ",maxerr);
					//mexPrintf("%d ",j);
                }
            }
        }
		//mexPrintf("\n");
        errs[i] = maxerr;
        err = err+maxerr;
        flags[index[i]] = 1;
        if(maxerr<*threshold)
        {
//             *len = i;
             break;
        }
//        pser[i] = (1-err)/((1-lambda*(i+1)/M)*(1-lambda*(i+1)/M));
    }
     
//    tk = minind(pser,int(M/3));
//    maxerr = 1 - pser[tk]*((1-lambda*(tk+1)/M)*(1-lambda*(tk+1)/M));
//    for(i=0;i<tk+1;i++)  index[i] = nind[index[i]];
//    *len = tk;
	if(i<N)
		*len = i;
	else
		*len = --i;
    return maxerr;
};
    

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, mxArray *prhs[])
{
    double * Y,*threshold;
    double * P;
    size_t M,N;
    M = mxGetM(prhs[1]);
    N = mxGetN(prhs[1]);
    Y = mxGetPr(prhs[0]);
    P = mxGetPr(prhs[1]);
    threshold = mxGetPr(prhs[2]);
    double err0=0;
    int len=0;
    double * errs = (double*)mxCalloc(N,sizeof(double));
    int *index = (int*)mxCalloc(N,sizeof(int));  
//     double *pser = (double*)mxCalloc(N,sizeof(double));
//     err0=err_cal(Y, P, M, N, 0, index, &len, errs, pser);
    err0=err_cal(Y, P, M, N, 0, index, &len, errs,threshold);
    plhs[0] = mxCreateDoubleMatrix(1,len+1,mxREAL); 
    mxSetPr(plhs[0],errs);
    plhs[1] = mxCreateNumericMatrix(1,len+1,mxINT32_CLASS,mxREAL); 
    mxSetData(plhs[1],index); 
    return;
}

            
            
        
        
    