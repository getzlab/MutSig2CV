function ok = demand_file(fname)

if ~iscell(fname), fname = {fname}; end

ok = true(size(fname));
for i=1:numel(fname), if ~mod(i,100), fprintf('%d/%d ',i,numel(fname)); end
  if ~exist(fname{i},'file'), ok(i)=false; end
end, if numel(fname)>=100, fprintf('\n'); end
nnf = sum(~ok);
nf = sum(ok);

if nargout==0
  if nf>=10, fprintf('%d files found\n',nf); end
end

if nnf>0
  blank = 0;
  for i=1:numel(fname)
    if ~ok(i)
      if strcmp(fname{i},'')
        blank=blank+1;
      else
        if nnf<20
          fprintf('\tNot found: %s\n',fname{i});
        end
      end
    end
  end
  if blank>0, fprintf('\tNot found: <blank> (%d)\n',blank); end
  if nargout==0
    if nnf==1, error('1 file not found'); else error('%d files not found',nnf); end
  end
end

if nargout==0, clear ok; end
