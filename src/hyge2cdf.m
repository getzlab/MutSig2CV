function p=hyge2cdf(kv,n,k1,n1)

if numel(kv)>length(kv) || numel(n)>length(n) || numel(k1)>length(k1) || numel(n1)>length(n1)
  error('only works with scalars/vectors; not matrices');
end

p = nan(size(kv));

if length(n)==1 && length(k1)==1 && length(n1)==1

  % kv can be vector
  % other parameters are scalar constants

  for i=1:length(kv)
    p(i)=0;
    for k=0:kv(i)
      term=gammaln(n1+2)-gammaln(k1+1)-gammaln(n1-k1+1) + ...         %denominator beta function
           gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1) + ...         %binom. coeff
           gammaln(k1+k+1) + gammaln(n+n1-k-k1+1) - gammaln(n+n1+2);  %numerator beta function
      prob=exp(term);
      p(i)=p(i)+prob;
    end
  end

else

  % all parameters are vectors

  len = length(kv);
  if length(n)~=len || length(k1)~=len || length(n1)~=len, error('inconsistent parameter lengths'); end
  nv = n;
  k1v = k1;
  n1v = n1;

  for i=1:len
    p(i)=0;
    n = nv(i);
    k1 = k1v(i);
    n1 = n1v(i);
    for k=0:kv(i)
      term=gammaln(n1+2)-gammaln(k1+1)-gammaln(n1-k1+1) + ...
           gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1) + ...
           gammaln(k1+k+1) + gammaln(n+n1-k-k1+1) - gammaln(n+n1+2);
      prob=exp(term);
      p(i)=p(i)+prob;
    end
  end

end


p(p<0)=0;
p(p>1)=1;
