function S = orderfields_first(S,first_flds)
% orders the given fields first; appends the rest of the fields in existing order
%
% Mike Lawrence 2009-11-11

if ischar(first_flds), first_flds = {first_flds}; end

all_flds = fieldnames(S);

if ~isempty(setdiff(first_flds,all_flds)), error('Some of those fields don''t exist'); end

rest_flds = all_flds;
rest_flds(ismember(rest_flds,first_flds)) = [];

S = orderfields(S,[as_column(first_flds);as_column(rest_flds)]);

