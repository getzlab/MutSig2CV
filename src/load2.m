function X = load2(indir,flds)

if (nargin~=1 && nargin~=2) || nargout~=1 || ~ischar(indir), error('usage:  X = load2(indir)'); end

if ~exist(indir,'dir')
  error('directory %s does not exist',indir);
end

dall = direc(indir);

f1=[];
f1.file = direc([indir '/*.mat']);
f1 = parse_in(f1,'file',[indir '/(.*)\.mat$'],'field');
f1.load2_by_chunks = false(slength(f1),1);

f2=[];
f2.file = direc([indir '/*.save2_by_chunks']);
f2 = parse_in(f2,'file',[indir '/(.*)\.save2_by_chunks$'],'field');
f2.load2_by_chunks = true(slength(f2),1);

d = concat_structs({f1,f2});
d.datenum = get_filedatenum(d.file);
d = sort_struct(d,'datenum');

if isempty(d) && ~isempty(dall)
  error('Could not load properly: is this a loadM situation?');
end

if isempty(d) && isempty(dall)
  fprintf('Warning: empty directory.');
end

if ~isempty(d) && any(~ismember(dall(~grepm('md5$', dall)),d.file))
  fprintf('Warning: some extraneous files in this directory.');
end

X=[];

nd = slength(d);
for i=1:nd
  if exist('flds','var') && ~ismember(d.field{i},flds), continue; end
  if nd>1, fprintf('%d/%d ',i,nd); end
  clear tmp;

  if ~d.load2_by_chunks(i)
    load(d.file{i},'tmp');
    if ~exist('tmp','var'), error('%s does not contain "tmp" variable',d.file{i}); end
  else % load2_by_chunks
    indir2 = d.file{i};
    d2 = direc([indir2 '/chunk*.mat']);
    if isempty(d2), error('problem with load2_by_chunks'); end
    tmp = cell(length(d2),1);
    for j=1:length(d2)
      clear tmp2;
      load(d2{j},'tmp2');
      if ~exist('tmp2','var'), error('%s does not contain "tmp2" variable',d2{j}); end
      tmp{j} = tmp2;
    end
    tmp = cat(1,tmp{:});
  end

  X = setfield(X,d.field{i},tmp);
end,fprintf('\n');
