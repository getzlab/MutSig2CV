function S = rmfield_if_exist(S,f)

if ~iscell(f), f = {f}; end

for i=1:length(f)
  if isfield(S,f{i}), S = rmfield(S,f{i}); end
end

