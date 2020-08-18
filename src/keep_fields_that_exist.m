function S2 = keep_fields_that_exists(S,flds)
% keep_fields_that_exists(S,flds)
%
% given struct <S> and cell-array-of-strings <flds>,
% returns struct <S2> which has only those fields specified.
%
% Mike Lawrence 2010-02-17
S2=[];
if isempty(S), return; end
for i=1:length(flds)
  if isfield(S,flds{i})
    f = getfield(S,flds{i});
    S2 = setfield(S2,flds{i},f);
  end
end
