function [ct,ua,ub] = xcount(a, b, sort_by_frequency, truncate_length)
% xcount(a,b,sort_by_frequency)
%
% Finds the unique pairs across two arrays (a and b) and counts how many of each pair there
%
% If sort_by_frequency is 1 or -1 (reverse), then table is sorted.
%
% if truncate_length is specified, trims strings to specified length (to improve visual fit)
%
% 2008-06-16 Mike Lawrence
% 2012-12-12 speedup for working with numbers
% 2016-12-01 modified to return struct if nargout==1

if ~ismember(nargout,[0 1 3]), error('invalid nargout'); end

if size(a,1)==1, a=a'; end
if size(b,1)==1, b=b'; end

if length(a)~=length(b), error('Arrays must be same length'); end

% find unique elements and make table

[ua uai uaj] = nanunique(a);
nua = length(ua);

[ub ubi ubj] = nanunique(b);
nub = length(ub);

% make sure ua and ub are cell arrays of strings
% (if they aren't, convert them)

if ~iscell(ua)
  tmp=ua;
  ua = cell(length(tmp),1);
  for i=1:length(tmp)
    if ischar(tmp(i))
      ua(i) = {tmp(i)};
    else
      ua(i) = {num2str(tmp(i))};
end,end,end

if ~iscell(ub)
  tmp=ub;
  ub = cell(length(tmp),1);
  for i=1:length(tmp)
    if ischar(tmp(i))
      ub(i) = {tmp(i)};
    else
      ub(i) = {num2str(tmp(i))};
end,end,end

% if truncate_length specified, trim values
if exist('truncate_length','var')
  for i=1:length(ua)
    ua{i} = ua{i}(1:min(truncate_length*2,length(ua{i})));
  end
  for i=1:length(ub)
    ub{i} = ub{i}(1:min(truncate_length,length(ub{i})));
  end
end

% compute 2D histogram
try
  ct = hist2d_fast(uaj,ubj,1,nua,1,nub);
catch
  ct = zeros(nua,nub);
  for aj=1:nua
    for bj=1:nub
      ct(aj,bj) = sum(uaj==aj & ubj==bj);
    end
  end
end

% compute totals
tota = sum(ct,2);
totb = sum(ct,1);

% sort if requested (rows by tota, cols by totb)
if exist('sort_by_frequency', 'var') && sort_by_frequency
  if sort_by_frequency == 1
    [tmp orda] = sort(tota);
    [tmp ordb] = sort(totb);
  elseif sort_by_frequency == -1
    [tmp orda] = sort(tota, 'descend');
    [tmp ordb] = sort(totb, 'descend');
  end
  ct = ct(orda,ordb);
  tota = tota(orda);
  totb = totb(ordb);
  ua = ua(orda);
  ub = ub(ordb);
end

% make cell table, including totals

tbl = cell(nua+2,nub+2);
for aj=1:nua
  for bj=1:nub
    tbl{aj+1,bj+1} = num2str(ct(aj,bj));
  end
end

for aj=1:nua
  tbl{aj+1,1} = ua{aj};
  tbl{aj+1,end} = num2str(tota(aj));
end

for bj=1:nub
  tbl{1,bj+1} = ub{bj};
  tbl{end,bj+1} = num2str(totb(bj));
end

tbl{end,1} = 'TOTAL';
tbl{1,end} = 'TOTAL';
tbl{1,1} = '';
tbl{end,end} = num2str(sum(tota));

% measure columns

L=zeros(nua+2,nub+2);
for aj=1:nua+2
  for bj=1:nub+2
    L(aj,bj)=length(tbl{aj,bj});
  end
end

W=zeros(nub+2,1);
for bj=1:nub+2
  W(bj) = max(L(:,bj));
end

if nargout==0
  % print table
  fprintf('\n');
  for aj=1:nua+2
    if aj==nua+2, fprintf('\n'); end
    fprintf('    ');
    for bj=1:nub+2
      if bj==2, fprintf('  '); end
      if bj==nub+2, fprintf('    '); end
      fprintf(['%' num2str(W(bj)) 's  '], tbl{aj,bj});
    end
    fprintf('\n');
  end
  fprintf('\n');

  clear ct ua ub;

elseif nargout==3
  % return ct ua ub

elseif nargout==1
  % return struct
  Q=[];
  Q.name = ua;
  for i=1:length(ub), Q.(genvarname(ub{i}))=ct(:,i); end

  ct = Q;
  clear ua ub;

else
  error('invalid nargout');
end
