function Y = move_field_to_before(X,fld1,fld2)
% Y = move_field_to_before(X,fld1,fld2)

demand_fields(X,{fld1,fld2});
f = fieldnames(X);

for i=1:length(f)
  if strcmp(f{i},fld1), continue; end
  if strcmp(f{i},fld2)
    Y.(fld1) = X.(fld1);
  end
  Y.(f{i}) = X.(f{i});
end
