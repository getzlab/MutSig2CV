function h = hist2d(y,x,ybins,xbins)

y = as_column(y);
x = as_column(x);
if length(y)~=length(x), error('y and x need to be same length'); end
h = zeros(length(ybins),length(xbins));
for i=1:length(y)
  yidx = find(y(i)>=ybins,1,'last');
  if isempty(yidx), yidx=length(ybins); end
  xidx = find(x(i)>=xbins,1,'last');
  if isempty(xidx), xidx=length(xbins); end
  h(yidx,xidx)=h(yidx,xidx)+1;
end
