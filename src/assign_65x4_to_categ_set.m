function c = assign_65x4_to_categ_set(K)
% c = assign_65x4_to_categ_set(K)
%
% Given struct K specifying k categories as either:
%
%   k rows and the following fields:
%   left    = subset of 'ACGT', representing 5' base
%   right   = subset of 'ACGT', representing 3' base
%   from    = subset of 'AC', representing mutated base (after strand collapse)
%   change  = subset of 'tfs', representing Transition, Flip transversion, Skew transversion
%   type   = either "point" or "non-point"
%
%   LEGACY CASE:  if name or type contains "null" or "indel" (or "double_null"), then type is taken to be "non-point"
%
%  or
%   192 or 96 rows and the following fields:
%   trinuc  = trimer of 'ACGT' (if middle base is only C/A then assumes strand collapse)
%   newbase = one of 'ACGT'
%   categ   = 1-k
%
% Maps the set of 65 territory categories and 4 newbases
%    onto this reduced set of categories.
%
% Returns matrix c, with 65 rows, k columns, and 4 pages (one for each newbase);
%   cells with "1" are counted toward that category.
%   For indel/null categories, all pages (newbases) are set to the same value.
%
% Mike Lawrence 2010-01-27
%  --> modified 2014-10-22 to allow piecewise definition of categories


X = generate_categ_context65_names();
require_fields(X,{'num','name'});
X = sort_struct(X,'num');  % (already sorted, but doesn't hurt to make sure)
X = parse(X.name,'(.) in (.)_(.)',{'from','left','right'});

base = 'ACGT';
complement(base) = 'TGCA';

if isfield(K,'trinuc') && isfield(K,'newbase') && isfield(K,'categ') && (slength(K)==96 || slength(K)==192)
  % NEW STYLE
  K = make_numeric(K,'categ');
  if any(K.categ<1 | isnan(K.categ)), error('bad entries in "categ" field of K'); end
  if slength(K)==96
    K2 = K;
    K2.trinuc = rc(K2.trinuc);
    K2.newbase = rc(K2.newbase);
    K = concat_structs({K,K2});
  end
  if length(unique_combos(K.trinuc,K.newbase))~=192, error('problem with categs_file'); end

  X.trinuc = stringsplice([X.left X.from X.right]);
  K.xi = listmap(K.trinuc,X.trinuc);
  K.ni = listmap(K.newbase,{'A','C','G','T'});
  if any(isnan(K.xi)) || any(isnan(K.ni)), error('problem with K'); end

  nk = max(K.categ);
  c = zeros(65,nk,4);
  for i=1:192, c(K.xi(i),K.categ(i),K.ni(i))=1; end

else
  % OLD STYLE
  require_fields(K,{'left','right','from','change','type'});
  nk = slength(K);
  % legacy cases
  idx = grepi('indel|null',K.type,1);
  if ~isempty(idx), K.type(idx) = repmat({'non-point'},length(idx),1); end
  if isfield(K,'name')
    idx = grepi('indel|null',K.name,1);
    if ~isempty(idx), K.type(idx) = repmat({'non-point'},length(idx),1); end
  end
  
  whatchange('A','ACGT') = 'nstf';
  whatchange('C','ACGT') = 'snft';
  whatchange('G','ACGT') = 'tfns';
  whatchange('T','ACGT') = 'ftsn';
  
  c = zeros(65,nk,4);
  
  for x=1:65   % for each context65 
    from = X.from{x};
    if isempty(from), continue; end   % "N"
    left = X.left{x};
    right = X.right{x};
    if ismember(from,'GT')
      from = complement(from);
      left = complement(X.right{x});
      right = complement(X.left{x});
    end
    for k=1:nk
      % decide if the context belongs to this category
      if ismember(from,K.from{k}) && ismember(left,K.left{k}) && ismember(right,K.right{k})
        % yes, it does
        if strcmpi(K.type{k},'non-point')
          c(x,k,:) = 1;   % for non-point mutations (indel/null), "newbase" doesn't matter
        elseif strcmpi(K.type{k},'point')
          for n=1:4     % which newbase
            oldbase = X.from{x};
            newbase = base(n);
            change = whatchange(oldbase,newbase);
            if ismember(change,K.change{k})
              c(x,k,n) = 1;
            end
          end
        else
          error('Unknown type: %s',K.type{k});
        end
      end     
    end  % next k
  end  % next x

end


