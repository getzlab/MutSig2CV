function y = my_hygepdf(x,m,k,n)
%HYGEPDF Hypergeometric probability density function.
%   Y = HYGEPDF(X,M,K,N) returns the hypergeometric probability
%   density function at X with integer parameters M, K, and N.
%   Note: The density function is zero unless X is an integer.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   See also HYGECDF, HYGEINV, HYGERND, HYGESTAT, PDF.

%   Reference:
%      [1]  Mood, Alexander M., Graybill, Franklin A. and Boes, Duane C.,
%      "Introduction to the Theory of Statistics, Third Edition", McGraw Hill
%      1974 p. 91.

%   Copyright 1993-2005 The MathWorks, Inc.
%   $Revision: 2.11.2.5 $  $Date: 2005/11/18 14:28:05 $

% (removed all error checking--> this gives a 25% speedup and allows 2D- or 3D-vectorized use)

kx = my_gammaln(k+1)-my_gammaln(x+1)-my_gammaln(k-x+1);
mn = my_gammaln(m+1)-my_gammaln(n+1)-my_gammaln(m-n+1);
mknx = my_gammaln(m-k+1)-my_gammaln(n-x+1)-my_gammaln(m-k-(n-x)+1);
y = exp(kx + mknx - mn);

