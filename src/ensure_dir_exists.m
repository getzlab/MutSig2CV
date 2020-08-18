function ensure_dir_exists(dirname)

if ~iscell(dirname)
  if ~exist(dirname,'dir'), mkdir(dirname); end
else
  for i=1:numel(dirname)
    if ~exist(dirname{i},'dir'), mkdir(dirname{i}); end
  end
end

