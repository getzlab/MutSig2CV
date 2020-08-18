function m = listmap(a,b,rowsflag)

if ~exist('rowsflag','var'), rowsflag=false; end

if rowsflag
  if ~(isnumeric(a) && isnumeric(b) && ndims(a)==2 && ndims(b)==2 && size(a,2)==size(b,2))
    error('not an appropriate input for "rows" operation');
  end

  na = size(a,1);
  nb = size(b,1);
  
  long = (na>=1e5 && nb>=1e5);

  if long, fprintf(' [init]'); end
  x = [];
  x.val = [a;b];
  x.a = [true(na,1);false(nb,1)];
  x.idx = as_column([1:na 1:nb]);
  if long, fprintf(' [unique]'); end
  [tmp ui x.uj] = unique(x.val,'rows');
  if long, fprintf(' [sort]'); end
  x = sort_struct(x,{'uj','a'});
  if long, fprintf(' [loop]'); end
  m = nan(na,1);
  aidx = nan;
  bidx = nan;
  for i=1:slength(x)
    if i>1 && x.uj(i)~=x.uj(i-1), aidx=nan; bidx=nan; end
    if x.a(i), aidx=x.idx(i); else bidx=x.idx(i); end
    if ~isnan(aidx), m(aidx)=bidx; end
  end
  if long, fprintf(' [done]\n'); end
  return
end
    

na = length(a);
nb = length(b);

if ischar(a) && ischar(b) && nb<na
  % special case of mapping characters to a string
  % typically something like listmap('ACAGCTGACACTGT','ACGT')
  % (later methods haven't been tested on this)
  method = 1;
else
  if ischar(a), a={a}; na=1; end
  if ischar(b), b={b}; nb=1; end
  method = 4;    % much faster for huge cases
end

long = (na>=1e5 && nb>=1e5);

if method==4
  if long, fprintf(' [init]'); end
  x = [];
  x.val = [as_column(a);as_column(b)];
  x.a = [true(na,1);false(nb,1)];
  x.idx = as_column([1:na 1:nb]);
  if long, fprintf(' [unique]'); end
  [tmp ui x.uj] = unique(x.val);
  if long, fprintf(' [sort]'); end
  x = sort_struct(x,{'uj','a'});
  if long, fprintf(' [loop]'); end
  m = nan(na,1);
  aidx = nan;
  bidx = nan;
  for i=1:slength(x)
    if i>1 && x.uj(i)~=x.uj(i-1), aidx=nan; bidx=nan; end
    if x.a(i), aidx=x.idx(i); else bidx=x.idx(i); end
    if ~isnan(aidx), m(aidx)=bidx; end
  end
  if long, fprintf(' [done]\n'); end

elseif method==3
  error('no such method');

elseif method==2  % was in use til replaced with method 4 in Mar2014

  if long, fprintf(' [init]'); end
  m = nan(na,1);
  if long, fprintf(' [unique]'); end
  [a ai aj] = unique(a);
  if long, fprintf(' [intersect]'); end
  [c ia ib] = intersect(a,b);
  if long, fprintf(' [loop]'); end
  for i=1:length(c), if ~mod(i,1e5), fprintf('%d/%d ',i,length(c)); end
    m(aj==ia(i)) = ib(i);
  end, if length(c)>=1e5, fprintf('\n'); end
  if long, fprintf(' [done]\n'); end

elseif method==1  % original method

  m = nan(na,1);
  [u ui uj] = unique(a);
  if ischar(u)
    for i=1:length(u)
      idx = find(b==u(i));
      if isempty(idx), idx = NaN; end
      m(uj==i) = idx;
    end
  else
    for i=1:length(u)
      idx = find(strcmp(u{i},b),1);
      if isempty(idx), idx = NaN; end
      m(uj==i) = idx;
    end
  end  
end
