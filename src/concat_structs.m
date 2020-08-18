function S = concat_structs(slist)
%
% concat_structs(slist)
%
% slist should be a cell array of structs,
%   all of which have the same fields
% 
% returns a new struct in which each field is a concatenation of
%   the corresponding fields from each of the input structs
%
% Mike Lawrence 2008-2010
%

if ~iscell(slist)
  error('input should be a cell array of structs');
end

if length(slist)==0
  S = {};
elseif length(slist)==1
  S = slist{1};
else

  ns = numel(slist);
  allflds = cell(ns,1);
  for i=1:ns
    if isempty(slist{i}), continue; end
    if ~isstruct(slist{i}), error('input should be a cell array of structs'); end
    allflds{i} = fieldnames(slist{i});
  end
  allflds = cat(1,allflds{:});
  [flds ai aj] = unique_keepord(allflds);
  h = histc(aj,1:length(flds));
  if ~all(h==ns)
    count(allflds);
    error('all input structs must have same fields.  use concat_structs_keep_all_fields.');
  end
  [tmp ord] = sort(ai);

  rflag = false;
  S = [];
  for fno=1:length(flds)
    type = nan(ns,1);
    f = cell(ns,1);
    for i=1:ns
      if isempty(slist{i}), continue; end
      f{i} = getfield(slist{i},flds{fno});
      if isnumeric(f{i}) || islogical(f{i}), type(i) = 1;
      elseif iscell(f{i}), type(i) = 2;
      else type(i) = 3;
      end
    end
    if any(type==3), error('Incompatible type encountered in %s',flds{fno}); end
    if any(type==1) && any(type==2)    
       if ~rflag, fprintf('Reconciling mixed cell+numeric fields for %s\n',flds{fno}); rflag=true; end
       for i=1:ns
         if type(i)==2
           try
             f{i} = str2double(f{i});
           catch me
             error('Unable to resolve mixed cell+numeric case encountered in %s',flds{fno});
    end,end,end,end
    S = setfield(S, flds{fno}, cat(1,f{:}));
  end

end
