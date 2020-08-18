/* projection_1d_convolutions_fast.c

   sped up version of projection_1d_convolutions.m

   function p = projection_1d_convolutions(Sdeg,Pdeg,score_obs,numbins,H,newH)

   input:
   Sdeg is score of each 1D degree (np,ncat+1):  MUST BE INTEGERS
   Pdeg is probability of each 1D degree (np,ncat+1)
   score_obs is observed score to calculate p-value for
   numbins is typically score_obs + ~10
   H is an externally allocated matrix (numbins,1)
   newH is an externally allocated matrix (numbins,ncols)

   returns:
   pmax = upper bound on p-value, i.e. the standard p-value usually calculated
   pmin = lower bound on p-value, going "one discrete step further"
 
   To compile: (from Matlab prompt)
   >> cd /xchip/cga2/lawrence/cga/trunk/matlab/seq
   >> mex projection_1d_convolutions_fast.c

*/

#include <string.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{
  double *Sdeg,*Pdeg;
  int *size;
  double *H, *newH;
  double score_obs;
  long numbins;
  long np,ncat,ncols,nnn;
  long i,j,k,p,d,idx,idx2;
  double odouble;
  long o;
  double pmax,pmin;
  double *returnval;

  if (nrhs!=6) mexErrMsgTxt("requires six input arguments: Sdeg,Pdeg,score_obs,numbins,H,newH");
  if (nlhs>2) mexErrMsgTxt("too many outputs requested");

  score_obs = mxGetScalar(prhs[2]);
  numbins = mxGetScalar(prhs[3]);

  size = (int *)mxGetDimensions(prhs[0]);
  np = size[0];
  ncat = size[1]-1;

  ncols = (ncat+1);

  /* binsize is fixed at 1 */

  size = (int *)mxGetDimensions(prhs[1]);
  if (size[0]!=np || size[1]-1!=ncat)
    mexErrMsgTxt("Pdeg should be same size as Sdeg");

  size = (int *)mxGetDimensions(prhs[4]);
  if (size[0]!=numbins) mexErrMsgTxt("H size should be numbins");
  
  size = (int *)mxGetDimensions(prhs[5]);
  if (size[0]!=numbins || size[1]!=ncols) mexErrMsgTxt("newH size should be (numbins,ncols)");

  Sdeg = mxGetPr(prhs[0]);
  Pdeg = mxGetPr(prhs[1]);
  H = mxGetPr(prhs[4]);
  newH = mxGetPr(prhs[5]);

  nnn = numbins*ncols;

  /* initial condition: all probability is in first bin */
  for (i=1;i<numbins;i++) H[i] = 0; 
  /* memset(H,0,8*numbins); */   /* using memset is ~3x slower! */
  H[0] = 1;

  /* sequential convolution */
  for (p=0;p<np;p++) {
    for (i=0;i<nnn;i++) newH[i] = 0; 
    /*memset(newH,0,8*nnn);*/
    idx = p;
    for (d=0;d<=ncat;d++) {
      odouble = Sdeg[idx];
      o = (odouble+0.5); /* round to nearest integer*/
      if (o<numbins) {
	idx2 = (d*numbins)+o;
	for (k=0; k<numbins-o; k++) newH[idx2+k] = Pdeg[idx] * H[k];
      } /* if */
      idx += np;
    } /* next d */
    /* sum newH across columns to get H */
    for (i=0;i<numbins;i++) {
      H[i] = 0;
      idx = i;
      for (j=0;j<ncols;j++) {
	H[i] += newH[idx];
	idx += numbins;
      }
    }
  } /* next patient */

  /* calculate p-value */


  /*
  pval = 1;
  for(i=0;i<numbins;i++) pval -= H[i];
  */

  pmax = 1;
  pmin = 1;

  for(i=0;i<numbins;i++) {
    if (H[i]>0) {
      pmin -= H[i];
      if (i+1>score_obs) break;
      pmax -= H[i];
    }
  }

  if (pmax<0) pmax=0;
  if (pmin<0) pmin=0;

  if (nlhs>=1) {
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    returnval = mxGetPr(plhs[0]);
    returnval[0] = pmax;
  }

  if (nlhs>=2) {
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    returnval = mxGetPr(plhs[1]);
    returnval[0] = pmin;
  }

}


