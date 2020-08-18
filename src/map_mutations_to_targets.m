function map = map_mutations_to_targets(x,T,P)

if ~exist('P','var'), P=[]; end

P = impose_default_value(P,'enforce_target_list_inclusive_of_noncoding',false);
P = impose_default_value(P,'enforce_target_list_footprint_table',[]);        % (only relevant if P.enforce_target_list_inclusive_of_noncoding==true)
P = impose_default_value(P,'enforce_target_list_footprint_flank_size',1e5);  % (only relevant if P.enforce_target_list_inclusive_of_noncoding==true)
% _____
% NOTE: if P.enforce_target_list_inclusive_of_noncoding==true
%       then the returned "map" values (pointers into T) serve only to indicate what gene the mutation is in.

fprintf('Mapping mutations to targets: ');

if P.enforce_target_list_inclusive_of_noncoding
  if isempty(P.enforce_target_list_footprint_table)
    error('cannot include noncoding mutations, because no footprint_table was provided');
  else
    F = P.enforce_target_list_footprint_table;
    if ischar(F), F=load_struct(F); end
    F = require_fields_with_convert(F,{'gene','chr','start','end'},{'name','chr','tx_start','tx_end'});
    F.chr = convert_chr(F.chr,P);
    F = make_numeric(F,{'start','end'});
    F.gene_idx = listmap(F.gene,T.gene);
    fprintf('\tIncluding noncoding mutations.\n');
  end
end

T.chr = convert_chr(T.chr,P);
T = make_numeric(T,{'start','end'});
x.chr = convert_chr(x.chr,P);
if ~isfield(x,'start') && isfield(x,'pos'), x.start=x.pos; end
if ~isfield(x,'end') && isfield(x,'pos'), x.end=x.pos; end
x = make_numeric(x,{'start','end'});

map = nan(slength(x),1);
for chr=1:24,fprintf('chr%d ',chr);
  tidx = find(T.chr==chr);
  xidx = find(x.chr==chr);
  for i=1:length(xidx), j = xidx(i);
    idx = tidx(find(x.start(j)<=T.end(tidx) & x.end(j)>=T.start(tidx),1));
    if ~isempty(idx), map(j) = idx; end
  end
  if P.enforce_target_list_inclusive_of_noncoding
    fidx = find(F.chr==chr);
    margins = [0 10.^(2:ceil(log10(P.enforce_target_list_footprint_flank_size)))]; margins(end)=P.enforce_target_list_footprint_flank_size;
    if isempty(margins), margins=[0]; end
    for margin = margins
      for i=1:length(fidx), j = fidx(i);
        map(xidx(isnan(map(xidx)) & x.start(xidx)<=F.end(j)+margin & x.end(xidx)>=F.start(j)-margin)) = F.gene_idx(j);
      end
    end
  end
end, fprintf('\n');

