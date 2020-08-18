function D = mapacross(A,B,C,filler)
% D = mapacross(A,B,C[,filler])
%
% for each item in A,
%   looks for that item in B.
%   if found, returns the corresponding item in C;
%   if not found, returns "filler".
% returns the items as D.
%
% Mike Lawrence 2009-07-14

if nargin~=3 && nargin~=4
  error('usage: D = mapacross(A,B,C[,filler])');
end

if ischar(A), single_string_flag=true; A={A}; else single_string_flag=false; end

if ~isnumeric(A) && ~iscell(A) && ~islogical(A)
  error('A should be numeric, cell, or logical');
end
if ~isnumeric(B) && ~iscell(B) && ~islogical(B)
  error('B should be numeric, cell, or logical');
end
if ~isnumeric(C) && ~iscell(C) && ~islogical(C)
  error('C should be numeric, cell, or logical');
end

if length(B) ~= length(C), error('length(B) ~= length(C)'); end

idx = listmap(A,B);
idx2 = find(~isnan(idx));

if ndims(C)>2 || (size(C,1)>1 && size(C,2)>1)
  flag=true;
  dsz = size(C); dsz(1) = length(A);
else
  flag=false;
  dsz = [length(A) 1];
end

if iscell(C)
  if ~exist('filler','var'), filler = ''; end
  D = repmat({filler},dsz);
elseif isnumeric(C)
  if ~exist('filler','var'), filler = nan; end
  if size(filler,2)==1 && dsz(2)>1, filler = repmat(filler,1,dsz(2)); end
  D = repmat(filler,dsz(1),1);
elseif islogical(C)
  if length(idx2)~=length(idx), error('Missing values when using logical output: behavior undefined'); end
  D = false(dsz);
else
  error('unknown output type');
end

if flag
  if ndims(C)>7, error('>7-D not supported'); end
  D(idx2,:,:,:,:,:,:) = C(idx(idx2),:,:,:,:,:,:);
else
  D(idx2) = C(idx(idx2));
end

if single_string_flag && iscellstr(D), D=decell(D); end

