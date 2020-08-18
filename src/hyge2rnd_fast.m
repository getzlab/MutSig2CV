function out = hyge2rnd_fast(N,x,X,mm,nn)

if ~exist('mm','var'), mm=1; end
if ~exist('nn','var'), nn=mm; end

out = nan(mm,nn);

for mi=1:mm, for ni=1:nn
    
    r = rand;

    n = 0;
    z = hyge2pdf(n,N,x,X);
    
    while z<r
      n=n+1;
      z=z+hyge2pdf(n,N,x,X);
    end

    out(mi,ni)=n;
end,end


