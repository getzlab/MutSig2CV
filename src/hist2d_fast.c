/*    hist2d_fast.c

      Mike Lawrence 2010-06-20

      To compile: (from Matlab prompt)
      >> cd /xchip/cga2/lawrence/cga/trunk/matlab/mike
      >> mex hist2d_fast.c

*/

#include <string.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{
  long xlen, ylen, len;
  long firstx, lastx, firsty, lasty, nxvals, nyvals;
  double *x, *y, *xvals, *yvals, *h;
  long i, xi, yi, idx;

  if (nrhs!=2 && nrhs!=6) mexErrMsgTxt("requires either two or six input arguments: x, y, [firstx, lastx, firsty, lasty]");
  if (nlhs>1) mexErrMsgTxt("too many outputs requested");

  xlen = mxGetN(prhs[0]); if (xlen==1) xlen = mxGetM(prhs[0]);
  ylen = mxGetN(prhs[1]); if (ylen==1) ylen = mxGetM(prhs[1]);
  if (xlen!=ylen) mexErrMsgTxt("error: x and y are of different lengths");
  len = xlen;
  x = mxGetPr(prhs[0]);
  y = mxGetPr(prhs[1]);

  if (nrhs==6) {
    /* user supplied ranges */
    firstx = (long)mxGetScalar(prhs[2]);
    lastx = (long)mxGetScalar(prhs[3]);
    firsty = (long)mxGetScalar(prhs[4]);
    lasty = (long)mxGetScalar(prhs[5]);
  } else {
    /* compute ranges as min and max of x and y */
    firstx = *x;
    lastx = *x;
    firsty = *y;
    lasty = *y;
    for (i=1;i<len;i++) {
      if ((*(x+i))<firstx) firstx = (*(x+i));
      if ((*(x+i))>lastx) lastx = (*(x+i));
      if ((*(y+i))<firsty) firsty = (*(y+i));
      if ((*(y+i))>lasty) lasty = (*(y+i));
    }
    printf("x range = %d to %d;  y range = %d to %d\n",firstx,lastx,firsty,lasty);
  }

  nxvals = lastx-firstx+1;
  nyvals = lasty-firsty+1;
  if (nxvals<1 || nyvals<1) mexErrMsgTxt("cannot have firstx>lastx or firsty>lasty");

  /*    x = rows        */
  /*    y = columns     */

  plhs[0] = mxCreateDoubleMatrix(nxvals,nyvals,mxREAL);
  h = mxGetPr(plhs[0]);
  for(i=0;i<len;i++) {
    xi = (long)(*(x+i))-firstx;
    yi = (long)(*(y+i))-firsty;
    if (xi<0 || xi>=nxvals || yi<0 || yi>nyvals) continue;  /* ignore out-of-range values */
    idx = (nxvals*yi+xi);
    (*(h+idx))++;
  }

}


