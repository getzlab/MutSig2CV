function sz = filedatenum(fname)

if ischar(fname)
  d = dir(fname);
  sz = nan(length(d),1);
  for i=1:length(d)
    sz(i) = d(i).datenum;
  end
elseif iscell(fname)
  sz = nan(length(fname),1);
  for i=1:length(fname)
    d = dir(fname{i});
    if ~isempty(d), sz(i) = d(1).datenum; end
  end
end

