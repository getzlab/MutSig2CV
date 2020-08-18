function B = stringsplice(A,dim,sep)
% stringsplice(A,dim,sep)
%
% given a cell matrix of strings <A>
% and a dimension <dim> (default=1),
% concatenates strings along the dimension and returns a cell array of strings.
%
% Mike Lawrence 2009-03-03

if nargin==2
  if ischar(dim)
    sep=dim;
    clear dim;
  end
end

if nargin==3
  if ischar(dim) && isnumeric(sep)
    tmp = dim;
    dim = sep;
    sep = tmp;
  elseif ischar(dim) || isnumeric(sep)
    error('unable to parse input parameters');
  end
end

if ~exist('dim','var'), dim=1; end
if ~exist('sep','var'), sep=''; end

if dim==1, B = cell(size(A,1),1);
else B = cell(1,size(A,2)); end

for i=1:length(B), if ~mod(i,1e5), fprintf('%d/%d ',i,length(B)); end
  if dim==1, B{i} = concat(A(i,:),sep);
  else B{i} = concat(A(:,i),sep); end
end, if i>=1e5, fprintf('\n'); end


