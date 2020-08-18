/*    count_overlaps_fast2

      input = column of "key"  (for example, can be chr*1e9+start)... must be SORTED by this value
              column of "pat_idx"
      output = patient-patient matrix, listing number of overlapping keys per patient pair

      To compile: (from Matlab prompt)
      >> cd /cga/tcga-gsc/home/lawrence/cgal/trunk/matlab/seq
      >> reuse .matlab-2013a
      >> mex <function_name>.c

*/

#include <string.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{
  double *key,*pat,*out;
  long npat,n1,n2,n,i,j,a,b,pa,pb,p;
  double keyi;

  if (nrhs!=2) mexErrMsgTxt("usage: count_overlaps_fast2(key,pat_idx)");

  /* input */
  key = (double *)mxGetData(prhs[0]);
  pat = (double *)mxGetData(prhs[1]);
  n1 = mxGetN(prhs[0]); if (n1==1) n1 = mxGetM(prhs[0]);
  n2 = mxGetN(prhs[1]); if (n2==1) n2 = mxGetM(prhs[1]);
  if (n1!=n2) mexErrMsgTxt("ERROR: first and second arguments must be same length");
  n = n1;

  /* find npat, and ensure key is sorted */
  npat = 1;
  for (i=0;i<n;i++) {
    if (i>0 && (*(key+i))<(*(key+i-1))) mexErrMsgTxt("key must be sorted");
    p = (*(pat+i));
    if (p<1) mexErrMsgTxt("pat_idx must be 1 or greater");
    if (p>npat) npat=p;
  }

  /* output */
  plhs[0] = mxCreateDoubleMatrix(npat,npat,mxREAL);
  out = mxGetPr(plhs[0]);

  /* count overlaps */
  i = 0;
  while(i<n) {
    j = i;
    keyi = (*(key+i));
    while (j+1<n && (*(key+j+1))==keyi) j++;
    if (j>i) {
      for (a=i;a<=j;a++) {
	pa = (*(pat+a))-1;
	for (b=i;b<=j;b++) {
	  pb = (*(pat+b))-1;
	  (*(out+pb*npat+pa))++;
	}
      }
    }
    i=j+1;
  }

  /* clear the diagonal */
  for (p=0;p<npat;p++) (*(out+p*npat+p)) = 0;

}


