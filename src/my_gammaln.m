function x = my_gammaln(a)

global cache

maxa = max(a(:));

if isempty(cache)
  cache = gammaln(1:maxa);
end

nc = length(cache);
if maxa>nc
  cache = [cache gammaln(nc+1:maxa)];
end

a = ceil(a);
a(a<1)=1;

x = nan(size(a));
x(:) = cache(a);

