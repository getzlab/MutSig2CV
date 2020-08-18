function M = prepare_for_calcsig_v1(M,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'use_bagel','true');
P = impose_default_value(P,'use_prior','true');

% data about category-specific mutation rate
% (exclude flanking mutations, because we don't yet have a reliable estimate of territory)
M.cat.ntot = histc(M.mut.categ(M.mut.is_coding),1:slength(M.cat));
M.cat.Ntot = sum((M.cov.gene_sil_repcov + M.cov.gene_non_repcov),1)' * slength(M.pat);
M.cat.mu = (M.cat.ntot./M.cat.Ntot);
x_tot = sum(M.cat.ntot);
X_tot = M.cat.Ntot(end);
mu_tot = x_tot/X_tot;
M.cat.murel = M.cat.mu/mu_tot;

% data about sample-specific mutation rate
% (exclude flanking mutations, because we don't yet have a reliable estimate of territory)
M.pat.ntot = histc(M.mut.pat_idx(M.mut.is_coding),1:slength(M.pat));
M.pat.Ntot = repmat(sum(M.cov.gene_sil_repcov(:,end) + M.cov.gene_non_repcov(:,end),1),1,slength(M.pat))';
M.pat.mu = (M.pat.ntot./M.pat.Ntot);
mu_tot = (sum(M.pat.ntot)/sum(M.pat.Ntot));
M.pat.murel = M.pat.mu/mu_tot;

% product of sample- and category-specific marginals
murel_catpat = bsxfun(@times,M.cat.murel',M.pat.murel);

% data about gene-specific mutation ratex_gene = M.gene.nsil + M.gene.nflank;
M.gene.x = M.gene.nsil + M.gene.nflank;
M.gene.X = M.gene.Nsil + M.gene.Nflank;
[M.gene.rate M.gene.ci_low M.gene.ci_high] = binofit_2d(M.gene.x,M.gene.X);
globalrate = sum(M.gene.x)/sum(M.gene.X);
M.gene.f = M.gene.ci_low/globalrate;

if P.use_bagel
  M.gene.x_bagel = M.gene.nfit - M.gene.x;
  M.gene.X_bagel = M.gene.Nfit - M.gene.X;

  % compute dependence of "w" upon "f"
  [f_sort ord] = sort(M.gene.f);
  binsize = 200;   % genes per bin
  nbins = ceil(slength(M.gene)/binsize);
  div = unique(f_sort(1:binsize:end));
  bin=[];bin.min = div; bin.max = [div(2:end);max(M.gene.f)+1];
  nbins = slength(bin);
  ww = 0:0.01:1;
  l=nan(length(ww),1);
  M.gene.w_bagel = nan(slength(M.gene),1);
  for i=1:nbins
    gidx = find(M.gene.f>=bin.min(i) & M.gene.f<bin.max(i));
    bin.gidx{i,1} = gidx;
    bin.ngenes(i,1) = length(gidx);
    for j=1:length(l)
      l(j)=bbinologlik(M.gene.x(gidx),M.gene.X(gidx),ww(j)*M.gene.x_bagel(gidx)+1,ww(j)*(M.gene.X_bagel(gidx)-M.gene.x_bagel(gidx))+1);
    end
    [tmp idx] = max(l);
    bin.w(i,1) = ww(idx);
    M.gene.w_bagel(gidx) = ww(idx);
  end

  % combine gene + bagel into a single x+X:
  M.gene.x_sphere = M.gene.x + M.gene.w_bagel.*M.gene.x_bagel;
  M.gene.X_sphere = M.gene.X + M.gene.w_bagel.*M.gene.X_bagel;

else % no bagel

  M.gene.x_sphere = M.gene.x;
  M.gene.X_sphere = M.gene.X;
  
end

% add prior from what the per-gene rates are like
if P.use_prior
  [M.prior.a,M.prior.b]=mle_beta(M.gene.x,M.gene.X);
  fprintf('Priors:  a = %f   b = %f\n',M.prior.a,M.prior.b);
  M.gene.x_sphere = M.gene.x_sphere + (M.prior.a-1);
  M.gene.X_sphere = M.gene.X_sphere + (M.prior.b+M.prior.a-2);
end

%% scale by mu_cs
M.gene.X_work = repmat(M.gene.X_sphere,[1 slength(M.cat) slength(M.pat)]);
M.gene.x_work = bsxfun(@times,M.gene.x_sphere,shiftdim(murel_catpat',-1));

% actual nonsilent mutations and territory
M.gene.N_work = repmat(M.cov.gene_non_repcov,[1 1 slength(M.pat)]);
M.gene.n_work = nan(slength(M.gene),slength(M.cat),slength(M.pat));
fprintf('patient:');
for i=1:slength(M.pat), if ~mod(i,10), fprintf(' %d/%d',i,slength(M.pat)); end
  midx = find(M.mut.pat_idx==i & M.mut.is_coding & ~M.mut.is_silent);
  M.gene.n_work(:,:,i) = hist2d_fast(M.mut.gene_idx(midx),M.mut.categ(midx),1,slength(M.gene),1,slength(M.cat));
end, fprintf('\n');

function M = calcsig_v1(M,P)

demand_field(M.gene,{'name','n_work','N_work'});

P.gene_names = M.gene.name;

M.gene.ntot = sum(M.gene.n_work,3);
M.gene.Ntot = sum(M.gene.N_work(:,end,:),3);

M.gene = rmfield_if_exist(M.gene,{'p','q','rank'});
M.gene = orderfields_first(M.gene,{'name'});

if isfield(M.gene,'X_work') && isfield(M.gene,'x_work') && ~isfield(M.gene,'mu_work')
  M.gene.p = calculate_significance(M.gene.N_work,M.gene.n_work,{M.gene.X_work,M.gene.x_work},P);
elseif ~isfield(M.gene,'X_work') && ~isfield(M.gene,'x_work') && isfield(M.gene,'mu_work')
  M.gene.p = calculate_significance(M.gene.N_work,M.gene.n_work,M.gene.mu_work,P);
else
  error('M.gene needs either mu_work OR X_work, x_work');
end

M.gene.q = calc_fdr_value(M.gene.p);

if isfield(M.gene,'effect')
  [tmp ord] = sort_struct(M.gene,{'p','effect'},[1 -1]);
else
  [tmp ord] = sort_struct(M.gene,'p');
end

[tmp M.gene.rank] = sort(ord);



