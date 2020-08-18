function out = collapse_192_to_categ_set(in,K,dim)
% c = collapse_192_to_categ_set(in,K,dim)
%
% Given 
% "in", a matrix with 192 rows
% (or 192 columns, with dim=2, etc...)
%
% and K, a set of k categories with the following fields:
%   left   = subset of 'ACGT', representing 5' base
%   right  = subset of 'ACGT', representing 3' base
%   base   = subset of 'ACGT', representing mutated base
%            NOTE: if from contains only A&C, then will assume we are doing strand collapse
%            NOTE: "base" can also be called "from"
%   change = subset of 'tfs', representing Transition, Flip transversion, Skew transversion
%     OR
%   newbase= subset of 'ACGT' but not <base>, representing the base changed to
%            if length(base)>1, then newbase cannot be used.
%
% and the set of 192 elementary mutation types
%   left   = one of 'ACGT', representing 5' base
%   right  = one of 'ACGT', representing 3' base
%   base   = one of 'ACGT', representing mutated base
%   newbase= one of 'ACGT' but not <base>, representing the base changed to
% NOTE: will assume the conventional default ordering.
%
% Maps in to out via mapping of C192 onto K.
%
% Returns vector out, with the 192 collapsed (by summing).
%
% NOTE: nulls/indels are not considered here.
% NOTE: gives an error if the categories in K are not surjective onto C192
%
% Mike Lawrence 2013-04-13

if nargin~=2 && nargin~=3, error('requires two inputs'); end

insz = size(in);
if exist('dim','var')
  if dim<1 || dim~=round(dim) || dim>length(insz), error('invalid dim'); end
else
  dim=1;
end
if insz(dim)~=192, error('require 192 elements in dimension %d',dim); end

map = assign_192_to_categ_set(K);
nk = max(map);
outsz = insz;
outsz(dim) = nk;
out = nan(outsz);

for k=1:nk
  j = (map==k);
  idxstr = repmat(':,',1,length(insz)); idxstr(end)=[];
  outidxstr = idxstr; outidxstr(2*dim-1)='k';
  inidxstr = idxstr; inidxstr(2*dim-1)='j';
  cmd = ['out(' outidxstr ') = sum(in(' inidxstr '),' num2str(dim) ');'];
  eval(cmd);
end


