function out = binornd_fast(N,p,mm,nn)

%if p>1, error('illegal p'); end
if ~exist('mm','var'), mm=1; end
if ~exist('nn','var'), nn=mm; end

out = nan(mm,nn);

for mi=1:mm, for ni=1:nn
    
    r = rand;

    n = 0;
    x = binopdf(n,N,p);
    
    while x<r
      %  fprintf('r=%f n=%d x=%f\n',r,n,x);
      n=n+1;
      x=x+binopdf(n,N,p);
    end

    out(mi,ni)=n;
end,end


