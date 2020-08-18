function require_fields(T,fields)

if ~iscell(fields)
   fields = {fields};
end

for i=1:length(fields)
  if ~isfield(T,fields{i})
    error(['Structure is missing required field "' fields{i} '"']);
  end
end

end
