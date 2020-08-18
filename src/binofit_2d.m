function [rate ci_low ci_high] = binofit_2d(n,N,varargin)

if ndims(n)~=ndims(N), error('n and N must have same dimensions'); end
if ~all(size(n)==size(N)), error('n and N must have same dimensions'); end
if ndims(n)>10, error('Does not support (N>10)-dimensional matrices'); end

rate = nan(size(n));
ci_low = nan(size(n));
ci_high = nan(size(n));

if any(N(:)<n(:))
  n = min(n,N);
%  fprintf('binofit: reducing n where n<N\n');
end

if any(n(:)~=round(n(:)))
  fprintf('binofit: rounding n\n');
  n = round(n);
end

if any(N(:)~=round(N(:)))
  fprintf('binofit: rounding N\n');
  N = round(N);
end

for i=1:size(n,2) 
  for j=1:size(n,3)
    for k=1:size(n,4)
      for l=1:size(n,5)
        for m=1:size(n,6)
          for o=1:size(n,7)
            for p=1:size(n,8)
              for q=1:size(n,9)
                for r=1:size(n,10)
                  [rate(:,i,j,k,l,m,o,p,q,r) c] = binofit(n(:,i,j,k,l,m,o,p,q,r),N(:,i,j,k,l,m,o,p,q,r),varargin{:});
                  ci_low(:,i,j,k,l,m,o,p,q,r) = c(:,1);
                  ci_high(:,i,j,k,l,m,o,p,q,r) = c(:,2);
end,end,end,end,end,end,end,end,end

if nargout==2
  ci_low = cat(ndims(n)+1,ci_low,ci_high);
  clear ci_high;
end





