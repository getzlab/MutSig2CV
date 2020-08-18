function S = rename_field(S, oldname, newname)
%
% rename_field(S, oldname, newname)
%
% renames a field in structure S
%
% Mike Lawrence 2008-05-15
%
% Modified 2008-11-18
% 1. accepts lists of field names
% Modified 2009-10-06
% 1. fixed bug when oldname{i}==newname{i}
% Modified 2011-04-20
% 1. don't reorder fields


if iscell(oldname) && iscell(newname)
  if ~iscell(newname) || length(oldname)~=length(newname), error('lists must be same length'); end
elseif ~iscell(oldname) && ~iscell(newname)
  oldname = {oldname};
  newname = {newname};
else
  error('improper parameters');
end

flds = fieldnames(S);

for i=1:length(oldname)
  f = getfield(S, oldname{i});
  S = setfield(S, newname{i}, f);
  if ~strcmp(oldname{i},newname{i})
    S = rmfield(S, oldname{i});
  end
  idx = find(strcmp(flds,oldname{i}));
  if length(idx)~=1, error('unexpected behavior'); end
  flds{idx} = newname{i};
end

S = order_fields_first(S,unique_keepord(flds));
