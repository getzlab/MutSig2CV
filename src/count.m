function [b,u] = count(a, sort_by_frequency, varargin)

if ~exist('sort_by_frequency', 'var'), sort_by_frequency = false; end

if nargin>=2 && length(sort_by_frequency)>1
  % probably meant to call xcount
  fprintf('Did you mean "xcount"?\n');
  xcount(a, sort_by_frequency, varargin{:})
  return
end

%% finds the unique elements of array and counts how many of each there

if size(a,1)==1
   a = a';
end

[u ui uj] = nanunique(a);
nu = length(u);

% make sure _u_ is cell array of strings

if ~iscell(a)
  tmp=u;
  u = cell(length(tmp),1);
  for i=1:length(tmp)
    if ischar(tmp(i))
      u(i) = {tmp(i)};
    else
      u(i) = {num2str(tmp(i))};
    end
  end
end

b = zeros(nu,1);
for j=1:nu
    b(j) = length(find(uj==j));
end

if sort_by_frequency == 1
  [b ord] = sort(b);
  u = u(ord);
elseif sort_by_frequency == -1
  [b ord] = sort(b, 'descend');
  u = u(ord);
end

bb = cell(nu,1);
for j=1:nu
    bb{j} = b(j);
end

u = [u; '----TOTAL'];
bb = [bb; {length(a(:))}];
nu=nu+1;

L=zeros(nu,1);
for i=1:nu
 L(i)=length(u{i});
end
maxL = max(L);

if nargout==0
  fprintf('\n');
  f = ['    %' num2str(maxL) 's: [%d]\n'];
  for i=1:nu
    fprintf(f, u{i}, bb{i});
  end
end

if nargout<2, clear u; else u=u(1:end-1); end
if nargout<1, clear b; end
