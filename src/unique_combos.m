function [u ui uj] = unique_combos(varargin)
% [u ui uj] = unique_combos(varargin)

% FAST WAY

z = nargin;

params={};
if length(varargin)>1 && ischar(varargin{end}) && (strcmpi(varargin{end},'last')|strcmpi(varargin{end},'first'))
  params{end+1} = varargin{end};
  varargin(end)=[];
  z=z-1;
end

n = length(varargin{1});

v = cell(z,1); vi = cell(z,1); vj = nan(n,z);

for i=1:z
  if length(varargin{i})~=n
    error('inconsistent argument lengths');
  end
  [v{i} vi{i} vj(:,i)] = nanunique(varargin{i},params{:});
end

[w ui uj] = unique(vj,'rows',params{:});
nw = size(w,1);
u = cell(nw,1);
for i=1:nw
  u{i} = cell(1,z);
  for j=1:z
    u{i}{j} = varargin{j}(ui(i));
  end
end

return

% SLOW WAY:

f = cell(nargin,1);

for i=1:nargin
  if isnumeric(varargin{i})
    f{i} = num2cellstr(varargin{i});
  elseif iscellstr(varargin{i})
    f{i} = varargin{i};
  else
    error('Unsure how to handle argument #%d',i);
  end
end

q = stringsplice(f,1,'###');

[u ui uj] = unique(q);




