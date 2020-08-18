function Y = nansub(X,idx,filler)
% nansub(X,idx[,filler])
%
% given a linear array X and a set of indices idx,
% returns a subset of X where each element of the subset
% is taken from X according to idx.
%
% NOTE: if X is a matrix, then nansub returns rows from the matrix
%
% if any elements of idx are nan (or nonpositive), then the function returns
% the element specified as filler.
%
% if filler is unspecified, then the default is
%    nan (for numeric or logical)
%    {}  (for cell)
%
% Mike Lawrence 2009-07-15

if length(size(X))==2 && size(X,1)==1 && size(X,2)>1
%   fprintf('note: converting first argument to column vector\n');
  X = X';
end

if iscellstr(X) && size(X,1)==1 && size(X,2)>1
  X=X';
end

if islogical(X)
  type = 0;
elseif isnumeric(X)
  type = 1;
elseif iscell(X)
  type = 2;
else
  error('Unsupported array type');
end

if ~exist('filler','var')
  if type==0
    filler = false;
  elseif type==1
    filler = nan;
  elseif type==2
    filler = {''};
  else
    error('Inconsistent behavior with "type"');
  end
end

if type==0
  if ~islogical(filler)
    error('Inappropriate filler for logical array');
  end
elseif type==1
  if ~isnumeric(filler)
    error('Inappropriate filler for numeric array');
  end
elseif type==2
  if ischar(filler)
    filler = {filler};
  end
  if ~iscell(filler)
    error('Inappropriate filler for cell array');
  end
else
  error('Inconsistent behavior with "type"');
end

sz = size(X); sz(1) = length(idx);
Y = repmat(filler,sz);
idx2 = find(~isnan(idx) & idx>=1 & idx<=length(X));
Y(idx2,:,:,:,:) = X(idx(idx2),:,:,:,:);

