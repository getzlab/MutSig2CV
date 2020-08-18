function T = require_fields_with_convert(T,fields,altnames)

if ~iscell(fields)
   fields = {fields};
end

for i=1:length(fields)
  if ~isfield(T,fields{i})
    replaced = false;
    if ~iscell(altnames{i})
      if isfield(T,altnames{i})
        T = rename_field(T,altnames{i},fields{i});
        replaced = true;
      end
    else
      for j=1:length(altnames{i})
        if isfield(T,altnames{i}{j})
          T = rename_field(T,altnames{i}{j},fields{i});
          replaced = true;
          break;
        end
      end
    end
    if ~replaced
      error(['Structure is missing required field "' fields{i} '"']);
    end
  end
end
