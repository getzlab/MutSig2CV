function output = pr(varargin)
% print rows
%
% usage:    pr(col1,col2,col3,...coln[,idx])
%        OR pr(X,[{'fld1','fld2','fld3',...,'fldn'}[,idx]])
%              where X is a struct of arrays
%
% add '#' as first argument to print record numbers
%
% Mike Lawrence 2011-2012

PAGEWIDTH = 300;

v = varargin;
if isempty(v), error('no input'); end

output = [];

if nargout>1, error('too many output variables'); end
toscreen = (nargout==0);
tovar = (nargout==1);

print_record_numbers = false;
if ischar(v{1}) && strcmp(v{1},'#')
  print_record_numbers = true;
  v = v(2:end);
end

if isstruct(v{1})
  s = v{1};
  v = v(2:end);
  flds = fieldnames(s);
  if ~isempty(v) && (ischar(v{1}) || iscellstr(v{1}))
    %% detect if second argument is a subset of fieldnames to print
    if ischar(v{1}), ff={v{1}}; else ff=v{1}; end
    if any(ismember(ff,flds))
      ff(~ismember(ff,flds)) = [];
      flds = ff;
      v = v(2:end);
    end
  end
  %% if an entries-to-print index is specified, then save time and memory by restricting to the printed rows NOW
  if ~isempty(v) && length(v)==1 && (isnumeric(v{1})||islogical(v{1}))
    orig_idx_field = 'tmpfld3524598677298456';
    s.(orig_idx_field) = (1:slength(s))';
    s = reorder_struct(s,v{1});
    v={};  % pretend there was never an index
    % (this fixes the problem with huge delays)
  end
  %% extract fields (TO DO: get rid of this copying, use some kind of pointers instead)
  x = cell(length(flds),1);
  for i=1:length(flds),x{i}=getfield(s,flds{i}); end
  totlen = slength(s);
else
  % it's not a struct: keep all fields that have the same length as the first field 
  len = nan(length(v),1);
  for i=1:length(v)
    sz = size(v{i});
    if length(sz)>2
      len(i) = inf;
      % it's multidimensional
    elseif length(sz)==2 && sz(1)>1
      len(i) = sz(1);
    else
      % it's a vector
      len(i) = length(v{i});
    end
  end
  idx = find(len~=len(1),1);
  if isempty(idx)
    % all fields are the same length:
    % find out if the last argument could possibly be meant as an index
    lf = v{end};
    if (length(v)>1 && length(lf)==numel(lf)) && ...
       (islogical(lf) ||...
       (isnumeric(lf) && min(lf)==1 && max(lf)==len(1) && length(unique(lf))==len(1)))
      %% then it's probably a sort index
      x = v(1:end-1);
      v = v(end);
    else
      x = v;
      v = [];
    end
  else
    x = v(1:idx-1);
    v = v(idx:end);
  end
  totlen = len(1);
end

clear idx

% detect if last argument is an index
if ~isempty(v)
  targ = v{1};
  if islogical(targ), targ=find(targ); end
  if isnumeric(targ) && numel(targ)==length(targ) && ~all(isnan(targ))
    targ(isnan(targ))=[]; % remove NaNs
    if all(targ>=1 & targ<=totlen)
      idx = targ;
      v = v(2:end); % remove it from the arg list
    end
  end
  if ischar(targ)
    targ = {targ};
  end
  if iscellstr(targ)
    for fi=1:length(x)
      if iscellstr(x{fi})
        if length(targ)==1
          tmp = find(strcmp(targ,x{fi}));
        else
          tmp = listmap(targ,x{fi});
          tmp(isnan(tmp)) = [];
        end
        if ~isempty(tmp)
          idx = tmp;
          v = v(2:end); % remove it from the arg list
          break
end,end,end,end,end

% there are some remaining unparsed arguments: choose an appropriate error mesage
if ~isempty(v)
  if length(v)==1 && ischar(v{1})
    error('failed to match %s',v{1});
  elseif length(v)==1 && iscellstr(v{1})
    expr = 'failed to match ';
    for i=1:min(length(v{1}),3)
      if i>1, expr = [expr ', ']; end
      expr = [expr v{1}{i}];
    end
    if length(v{1})>3, expr = [expr ', ...']; end
    error(expr);
  elseif length(v)==1 && isempty(v{1})
    fprintf('nothing to print\n'); if toscreen ,clear output; end, return;
  elseif length(v)==1 && isnumeric(v{1})
    error('index out of range');
  else
    error('unknown argument structure / failed to match index expression');
  end
end

%%%%%%%%%%%%%%%

if ~exist('idx','var')
  if length(x{1})==numel(x{1})
    idx = 1:length(x{1});
  else  % 2-D
    idx = 1:size(x{1},1);
  end
end
idx=as_column(idx);

if exist('flds','var')
  flds = as_column(flds);
  idx=[0;as_column(idx)];
  have_headers = true;
else
  flds=num2cellstr((1:size(x,1))');
  have_headers = false;
end

if print_record_numbers
  idx_to_print = (1:size(x{1},1))';   % 2018-07-24  FIXED BUG
%  idx_to_print = idx(idx>0);
  if exist('orig_idx_field','var'), idx_to_print = s.(orig_idx_field)(idx_to_print); end
  if size(x,2)==1, x = [{idx_to_print};x]; else x = [{idx_to_print},x]; end
  flds = [{'#'};flds];
end

% PRE-FORMATTING

% (1) find multidimensional fields and convert them to well-aligned strings
%     (except 2D fields with too many columns)
multi_flag = false(length(x),1);
for c=1:length(x)
%  if size(x{c},1)>1 && size(x{c},2)>1 && size(x{c},2)<=50 && ndims(x{c})<3 && (isnumeric(x{c}) || islogical(x{c}))
% 2018-07-24  FIXED BUG that was printing only the first column of multi-column fields, in cases of printing only a single row
  if size(x{c},2)>1 && size(x{c},2)<=50 && ndims(x{c})<3 && (isnumeric(x{c}) || islogical(x{c}))
    multi_flag(c) = true;
    if islogical(x{c}) || all(x{c}(:)==round(x{c}(:)))
      fmt = '%d ';
    elseif mean(x{c}(:)==round(x{c}(:)))>0.5
      %fmt = '%.1d ';
      fmt = '%0.2f ';
    elseif mean(x{c}(:)>1)>0.5
      fmt = '%3f ';
    elseif mean(x{c}(:)>0.5)>0.5
      fmt = '%.1f ';
    elseif mean(x{c}(:)>0.01)>0.4
      fmt = '%.2f ';
    else
      fmt = '%.1d ';
    end
    tmp = num2str(x{c},fmt);
    x{c} = cell(size(tmp,1),1);
    for i=1:size(tmp,1)
      x{c}{i} = ['  ' tmp(i,:) '  '];
    end
  end
end

% FORMATTING

z = cell(length(idx),length(x));
fmt = '%f';
for j=1:length(idx)
  i=idx(j);
  for c=1:length(x)
    if i==0
      % HEADER ROW
      z{j,c} = flds{c};
    else
      % NOT HEADER ROW
      if size(x{c},1)>1 && size(x{c},2)>1
        % multidimensional (or 2D with too many columns)
        if size(x{c},2)<=50 && ((isnumeric(x{c}) || islogical(x{c})) && ndims(x{c})<3)
          fprintf('WARNING: should not be encountered, because now handled in pre-formatting\n');
          z{j,c} = num2str(x{c}(i,:));
        end
        z{j,c} = '{...}';
        continue
      else
        if isnumeric(x{c}) || islogical(x{c})
          val = x{c}(i);
          if islogical(val) || val==round(val)
            fmt = '%d';
          elseif abs(val)<0.0001 || abs(val)>1e6
            fmt = '%.2d';
          elseif abs(val)>1
            fmt = '%.1f';
          elseif abs(val)>=0.1
            fmt = '%.2f';
          elseif abs(val)>=0.01
            fmt = '%.3f';
          elseif abs(val)>=0.001
            fmt = '%.4f';
          elseif abs(val)>=0.0001
            fmt = '%.5f'; 
          end
        elseif iscell(x{c})
          val = x{c}{i};
          if ischar(val) && size(val,1)==1
            % just print string directly
            z{j,c} = val;
            continue;
          else
            z{j,c} = '{...}'; continue;
            %        subfprintf('Don''t know how to handle value in field %s:\n',flds{c});
            %        whos val
            %        return
          end
        else
          z{j,c} = '{...}'; continue;
          %      subfprintf('Don''t know how to handle field %s:\n',flds{c});
          %      field = x{c};
          %      whos field
          %      return
        end
        z{j,c} = sprintf(fmt,val);
      end
    end
  end
end

% get column lengths
w = nan(length(x),1);
for c=1:length(x)
  w(c) = max(cellfun(@length,z(:,c)));
end

% for long columns, pad the header names with "___" to make it clearer
if have_headers
  for c=1:length(x)
    headlen = length(z{1,c});
    collen = w(c);
    if collen>headlen+6 || multi_flag(c)
      z{1,c} = [z{1,c} repmat('_',1,collen-headlen)];
    end
  end
end

% choose number of "pages"
pagew = PAGEWIDTH;
nsp = 1;

ncol = length(x);
totw = sum(w) + nsp*(ncol-1);
npages = ceil(totw/pagew);

lastcol=0;
for pi=1:npages
  pw = 0;
  firstcol = lastcol+1;
  lastcol = firstcol;
  if firstcol>1
    pw=pw+w(1)+nsp;
  end
  while lastcol<ncol && pw+w(lastcol+1)+nsp<=pagew
    pw=pw+w(lastcol+1)+nsp;
    lastcol=lastcol+1;
  end

  if firstcol>ncol, break; end

  for j=1:length(idx)
    if firstcol>1, cols = [1 firstcol:lastcol]; else cols = firstcol:lastcol; end
    for c=cols
      txt = z{j,c};
      if j==1 && have_headers, ansi_start_underline(); end
      subfwritepad(txt,w(c));
      if j==1 && have_headers, ansi_stop_underline(); end
      if c<length(x), subfwritepad('',nsp); end
    end
    subfprintf('\n');
  end

if pi<npages, subfprintf('\n'); end
end

if toscreen
  clear output
end


  function subfprintf(varargin)
    txt = sprintf(varargin{:});
    if tovar, output = [output txt]; end
    if toscreen, fprintf(txt); end
  end

  function subfwritepad(txt,width)
    ns = width-length(txt);
    if tovar, output = [output txt repmat(' ',1,ns)]; end
    if toscreen, fwrite(1,[txt repmat(' ',1,ns)]); end
  end


end % main function
