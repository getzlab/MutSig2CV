function [u ui uj] = nanunique(a,varargin)
% [u ui uj] = nanunique(a)
%
% Mike Lawrence 2010-08-20

[u ui uj] = unique(a,varargin{:});

if isnumeric(a)
  idx = find(isnan(u));
  if length(idx)>1
    if idx(end)~=length(u) || any(diff(idx)~=1), error('expected unique to place all NaNs at end'); end
    uj(isnan(a)) = idx(1);
    u(idx(2:end)) = [];
    ui(idx(2:end)) = [];
  end
end
  
