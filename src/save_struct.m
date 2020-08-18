function F = save_struct(S,filename,noheader,deprecated_formatflag)
%
% save_struct(S, filename)
%
% writes a tab-delimited table from the given struct.
% writes a header line based on the field names.
%
% fields are written in the order in which they appear in the struct.
%
% Mike Lawrence 2008-06-20
% modified 2009-07-15 to handle empty S
% modified 2011-05-16 to handle extremely long structs (build as chunks)
% modified 2012-01-25 to restore correct handling of empty S

if nargin==0, error('requires an input argument'); end
if nargin>4, error('too many input arguments'); end
if nargin==4, fprintf('USE OF DEPRECATED FORMATFLAG\n'); end

if nargout==0 && nargin<2, error('requires at least 2 arguments'); end
if nargout==1 && (nargin<2 || isempty(filename))
  return_string_only = true;
else
  return_string_only = false;
end
if nargout>1, error('only a single output possible'); end

if exist('noheader','var') && ~isempty(noheader) && (~islogical(noheader) || ~(noheader==false))
  if (islogical(noheader) && noheader==true) || strncmpi(noheader,'no_header',9) || strncmpi(noheader,'noheader',8)
    noheader = true;
  else
    error('third parameter should be "no_headers", true, false, empty, or nothing.');
  end
else
  noheader = false;
end

if ~return_string_only
  out = ensure_fopen(filename,'wt');
end


if ~isstruct(S)
  if isempty(S)
    S = struct;
  else
    error('S should be a struct');
  end
end

fld = fieldnames(S);
nf = length(fld);
slen = slength(S);

% see if struct is empty
if slen==0
  if return_string_only
    F = '';
  else
    for i=1:nf
      fprintf(out,fld{i});
      if i==nf
        fprintf(out,'\n');
      else
        fprintf(out,'\t');
      end
    end
    fclose(out);
  end
  return
end


% see if struct is to long too handle all at once

chunksize = round(1e7/nf);
if slen>chunksize

  if return_string_only
    error('struct is too large to convert all at once in memory');
  end

  for st=1:chunksize:slen
    fprintf('HUGE STRUCT: SAVING CHUNK %d/%d\n', ceil(st/chunksize), ceil(slen/chunksize));
    en = min(slen,st+chunksize-1);
    Si = reorder_struct(S,st:en);
    if exist('deprecated_formatflag','var')
      F = save_struct(Si,[],noheader,deprecated_formatflag);
    else
      F = save_struct(Si,[],noheader);
    end
    fwrite(out,F);
    clear F;
    noheader = 'no_header'; % subsequent chunks omit header
  end

else   % struct is not too big to save all at once

  F = cell(1,nf*2);
  nr = -1;
  tt0 = tic;

  for f=1:nf
    tt1 = toc(tt0);
    if tt1>10
      if ~exist('flag01','var')
        fprintf('  [save_struct] ');
        flag01 = true;
      end
      fprintf('%d/%d ',f,nf);
    end
    C = getfield(S,fld{f});    % get column
                               % check for legal type
    if isempty(C)
      C = {};
    elseif isnumeric(C) || islogical(C)
      if ndims(C)>2
        fprintf('Field %s is multidimensional: skipping\n', fld{f});
        C = {};
      elseif size(C,2)>1
        fprintf('Field %s is not a column vector: skipping\n', fld{f});
        C = {};
      else
        if exist('deprecated_formatflag','var') && deprecated_formatflag     % convert to cell array of strings
          C = cellstr(num2str(C));
        else                       % note: "-" is important to avoid extra spaces in output
          if any(mod(C,1))
            C = cellstr(num2str(C,'%-d'));     % not all integers
          else
            C = cellstr(num2str(C,'%-.0f'));    % all integers
          end
        end
      end
    else  % column is a cell
      idx = find(cellfun(@iscell,C));
      if ~isempty(idx)
        fprintf('WARNING: Field %s contains %d entries that are cells:\n', fld{f}, length(idx));
        fprintf('         replacing these with "?????{cell}?????"\n');
        C(idx) = repmat({'?????{cell}?????'},length(idx),1);
      end
      idx = find(~cellfun(@isempty,C) & ~cellfun(@ischar,C));
      if ~isempty(idx)
        fprintf('WARNING: Field %s contains %d entries that are not chars:\n', fld{f}, length(idx));
        fprintf('         replacing these with "?????{unknown_format}?????"\n');
        C(idx) = repmat({'?????{unknown_format}?????'},length(idx),1);
      end
    end
    if isempty(C)
      if nr==-1, error('Problematic leftmost column'); end
      C = repmat({'...'},nr,1);
    end
    % check for consistent column length
    if nr==-1, nr=length(C); end
    if nr~=length(C), error('Field %s is a different length', fld{f}); end
    % add column title
    if ~noheader, C = [fld{f}; C]; end
    % add column to file
    F{f*2-1} = C;
    % add tab or newline
    F{f*2} = repmat({char(9+(f==nf))},length(F{f*2-1}),1);
  end

  if tt1>10, fprintf(' [collapse]'); end
  F = strcat(F{:});              % collapse columns to lines
  F = [F{:}];                    % collapse lines to file
  
  if ~return_string_only
    if tt1>10, fprintf(' [write]'); end
    fwrite(out,F);
  end

  if tt1>10, fprintf('\n'); end
end

if ~return_string_only
  fclose(out);
  clear F
end





