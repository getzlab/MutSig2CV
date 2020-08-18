function X = new_find_duplicate_samples(m,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'duplicate_samples_criterion',1);

if ~isfield(m,'pos') && isfield(m,'start'), m.pos = m.start; end
demand_fields(m,{'patient','chr','pos'});
has_type = isfield(m,'type');
has_clas = isfield(m,'classification');
if has_type
  fprintf('Comparing on the basis of coding mutations only...\n');
  m = reorder_struct_exclude(m,grepmi('intron|igr|utr|rna|flank',m.type));
end

% assign priority for choosing which mutation to keep, if there are duplicate patient-chr-start combinations
if has_type & has_clas
  m.is_ind = grepmi('DEL|INS',m.classification);
  m.is_mis = grepmi('missense',m.type);
  m.is_trunc = grepmi('nonsense|non.?stop|splice|read.?thr',m.type);
  m.is_syn = grepmi('silent|synon',m.type);
  m.priority = -m.is_syn + 0 + 1*m.is_mis + 2*m.is_trunc + 3*m.is_ind;
  m = sort_struct(m,'priority');
  m = rmfield(m,{'is_ind','is_mis','is_trunc','is_syn','priority'});
end

% keep at most one mutation (of highest priority) per patient-chr-start
[u ui uj] = unique_combos(m.patient,m.chr,m.pos);
m = reorder_struct(m,ui);
M = maf2M(m);
M.mut.chr = convert_chr(M.mut.chr);
M.mut = make_numeric(M.mut,'pos');
M.mut = reorder_struct_exclude(M.mut,isnan(M.mut.chr)|isnan(M.mut.pos));
M.mut.key = 1e9*M.mut.chr+M.mut.pos;
M.mut = sort_struct(M.mut,'key');

% find duplicate patients by mutational overlap
ni = count_overlaps_fast2(M.mut.key,M.mut.pat_idx);
ti = bsxfun(@max,M.pat.nmut,M.pat.nmut');
fi = ni./ti;

if P.duplicate_samples_criterion==1
  ov = (ni>=10 & fi>=0.1) | (ni>=3 & fi>=0.3);
elseif P.duplicate_samples_criterion==2
  ov = (ni>=10 & fi>=0.1) | (fi>=0.3);
else
  error('invalid P.duplicate_samples_criterion');
end

fprintf('%d patients involved in an overlap.\n',sum(sum(ov,1)>0));

% make list of mutationally overlapping patients
W = [];
[a b] = find(ov);
[idx] = find(ov);
W.pat1 = M.pat.name(a);
W.pat2 = M.pat.name(b);
W.ttype1 = M.pat.ttype(a);
W.ttype2 = M.pat.ttype(b);
W.nmut1 = M.pat.nmut(a);
W.nmut2 = M.pat.nmut(b);
W.ni = ni(idx);
W.fi = fi(idx);
W = reorder_struct(W,b>a);

W = sort_struct(W,'ni',-1);

% find all cliques of overlaps
q = scomponents(sparse(a,b,1));

U=[];
[U.num ui uj] = unique([a;b]);
U.name = M.pat.name(U.num);
U.ttype = M.pat.ttype(U.num);
U.nmut = M.pat.nmut(U.num);
U.clique = q(U.num);
[tmp tmp U.clique] = unique(U.clique);
U = sort_struct(U,'clique');
U = order_field_first(U,'clique');

if slength(U)>0
  fprintf('%d cliques of overlapping samples.\n',max(U.clique));
  fprintf('%d unique samples.\n',slength(M.pat)-slength(U)+max(U.clique));
end

% Choose patients to remove (in each clique, keep the patient with the most mutations)
z = sort_struct(U,'nmut');
[q qi qj] = unique(z.clique);
keep = z.name(qi);
drop = setdiff(z.name,keep);

%%%%% 
% report:

X=[];
X.pat = M.pat;
X.ni = ni;
X.ti = ti;
X.fi = fi;
X.ov = ov;
X.W = W;
X.U = U;
X.drop = drop;









