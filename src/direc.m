function t = direc(dirname)

if exist(dirname,'dir')
  if dirname(end)~='/', dirname = [dirname '/']; end
end

t = {};

if contains(dirname,'?') error('"?" wildcard not supported.  Use "*".'); end

tmp = rdir(dirname);
if isempty(tmp), return; end
t = list2cell(tmp.name);
t = remove_ansi(t);
t = grepv('^\.',t);



return

%% (below = obsolete)

if ~isempty(grep('\*.*\/.+',dirname))

  % wildcards in path structure:
  % use slow method
  fprintf('Getting directory via system call...\n');
  [status result] = system(['ls -1 ' dirname]);
  if status~=0
    fprintf('ERROR in direc: please help!\n');
    keyboard
  end
  ra = remove_ansi(result);
  t = text_to_lines(ra);

else

  % no wildcards in path structure:
  % use fast method
  t = {};
  h = parse(dirname,'^(.*/)?([^/]*)$',{'dir','mask'});
  tmp = dir(dirname);
  if isempty(tmp), return; end
  t = list2cell(tmp.name);
  t = remove_ansi(t);
  t = grepv('^\.',t);
  t = stringsplice([repmat(h.dir,length(t),1) t]);

end
