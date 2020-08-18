function n = replace_with_simulated_data(N,mu,P)

if ~exist('P','var'), P=[]; end

[ng ncat np] = size(N);

% is mu supplied as {X,x} for beta-binomial?
if isnumeric(mu)
  using_betabinom = false;
  if ndims(mu)~=ndims(N) || any(size(mu)~=size(N)), error('size mismatch'); end
elseif iscell(mu) && length(mu)==2
  using_betabinom = true;
  X = mu{1};
  x = mu{2};
  clear mu;
  if ndims(x)~=ndims(N) || any(size(x)~=size(N)), error('size mismatch'); end
  if ndims(X)~=ndims(N) || any(size(X)~=size(N)), error('size mismatch'); end
%  a = x+1;
%  b = X-x+1;
%  clear x X;
else
  error('unknown format for mu');
end

%%%%%%%%%% SIMULATOR

fprintf('SIMULATING DATA: ');

if isfield(P,'genes_to_analyze')
  gidx = listmap(P.genes_to_analyze,P.gene_names);
else
  gidx=1:ng;
end

n = nan(size(N));
for g=gidx, if ~mod(g,1000), fprintf('%d/%d ',g,ng); end
disp(g)
  for c=1:ncat
    for p=1:np
      if using_betabinom
%%        n(g,c,p) = bbinornd_fast(N(g,c,p),a(g,c,p),b(g,c,p));
       n(g,c,p) = hyge2rnd_fast(N(g,c,p),x(g,c,p),X(g,c,p));
      else
        n(g,c,p) = binornd_fast(N(g,c,p),mu(g,c,p));
      end
end,end,end, fprintf('\n');

keyboard


