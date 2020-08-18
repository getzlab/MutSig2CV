function buildno = interpret_build(build)

buildno = nan;

if isnumeric(build)
  if numel(build)>1, error('interpret_build not compatible with vector'); end
  if build==18 || build==36
    buildno = 18;
  elseif build==19 || build==37
    buildno = 19;
  elseif build==38
    buildno = 38;
  end
elseif ischar(build)
  if strcmp(build,'hg18') || strcmp(build,'36') || strcmp(build,'18')
    buildno = 18;
  elseif strcmp(build,'hg19') || strcmp(build,'37') || strcmp(build,'19')
    buildno = 19;
  elseif strcmp(build,'hg38') || strcmp(build,'38') || strcmpi(build,'GRCh38')
    buildno = 38;
  end
else
  fprintf('Unknown format of build\n');
  whos build;
end

if isnan(build)
  fprintf('Unknown build: ');
  disp(build);
end

