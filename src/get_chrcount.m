function ct = get_chrcount(build)
% ct = get_chrcount(build)

if ~exist('build','var'), build = []; end
if isempty(build)
  fprintf('get_chrcount: assuming human\n');
  ct = 24;
  return;
end

if isnumeric(build)
  if build==35, build=17; end
  if build==36, build=18; end
  if build==37, build=19; end
  build = ['hg' num2str(build)];
end

if strcmp(build,'hg19gencode'), build='hg19'; end

% see if "build" is actually a path to the annotation directory.
% separate into annotation directory and root if so
if build(1) == '/',
  tmp = regexp(build, '(/.*/)(hg\d+)/?$', 'tokens');
  dirstem = tmp{1}{1};
  build = tmp{1}{2};
else
  dirstem = '/xchip/cga/reference/annotation/db/ucsc/';
end
%tmp = parse(build,'^/.*/(hg\d+)$','build');
%if ~isempty(tmp.build{1})
%  build = tmp.build{1};
%end

% METHODS AVAILABLE:
% if ReferenceInfoObj has been initialized, use it.
% if ReferenceInfoObj has not been initialized, use legacy methods.

try

  maxnum = ReferenceInfoObj.getMaxNum(build);
  ct = 1;
  for i=1:maxnum
    if strcmp('D',ReferenceInfoObj.getUse(i,build)), ct=i; end
  end

catch me

  % legacy methods

  dr = [dirstem '/' build];

  fname = [dr '/cytoBand.txt'];
  if exist(fname,'file')
    C = load_struct_noheader(fname);
    
    if strcmp(C.col1{1},'chr'), C = reorder_struct(C,2:slength(C)); end
    
    uc=[]; uc.name = unique(C.col1);
    uc.name = regexprep(uc.name,'^chrM$','chr0');
    uc.name = regexprep(uc.name,'^chrMt$','chr0');
    uc.name = regexprep(uc.name,'^chrMT$','chr0');
    uc = parse_in(uc,'name','^chr(\d+)$','num',1);
    
    ct = slength(uc);
  else
    fname = [dr '/' build '_info.txt'];
    if exist(fname,'file')
      C = load_struct(fname);
      if isfield(C,'use'), C = reorder_struct(C,strcmp('D',C.use)); end
      C = make_numeric(C,'num');
      ct = max(C.num);
    else
      error('Don''t know how to find chrcount for build %s',build);
    end
  end

end






