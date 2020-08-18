function G = MutSig_2CV_v3_11_core(mutation_file, output_dir, mutation_blacklist_file, mutation_type_dictionary_file, coverage_models_mat_file, ...
                                      basewise_coverage_fwb_file, target_list_file, context_and_effect_fwb_file, context_and_effect_categs_file, ...
                                      covariates_file, conservation_fwb_file, FixedWidthBinary_jar_file, build, params_file)

if nargin<13, error('requires 13 input arguments'); end
if nargin>14, error('extraneous arguments'); end

% output version report textfile
MUTSIG_VERSION = '2CV v3.11'
outname = [output_dir '/MutSig_version.txt'];
save_textfile(['MutSig ' MUTSIG_VERSION], outname);

% ensure input files exist
demand_file(mutation_file);
if ~strcmp('',mutation_blacklist_file) & ~strcmp('none',mutation_blacklist_file), demand_file(mutation_blacklist_file); end
demand_file(mutation_type_dictionary_file);
demand_file(coverage_models_mat_file);
demand_file(basewise_coverage_fwb_file);
demand_file(target_list_file);
demand_file(context_and_effect_fwb_file);
demand_file(context_and_effect_categs_file);
demand_file(covariates_file);
demand_file(conservation_fwb_file);
demand_file(FixedWidthBinary_jar_file);

% add jar to java classpath
javaclasspath(FixedWidthBinary_jar_file);

% ensure output directory is ok
ede(output_dir);

P=[];
if exist('params_file','var') && ~strcmp(params_file,'') && ~strcmp(params_file,'none')
  P = process_params_file(P,params_file);
end

% MutSig parameters
P = impose_default_value(P,'enforce_target_list',true);
P = impose_default_value(P,'enforce_target_list_inclusive_of_noncoding',true);
P = impose_default_value(P,'enforce_target_list_footprint_flank_size', 3000);
P = impose_default_value(P,'infer_fields',true);
P = impose_default_value(P,'conform_to_maf_spec',false);
P = impose_default_value(P,'number_of_categories_to_discover',5);
P = impose_default_value(P,'gene_min_frac_coverage_required',0);
P = impose_default_value(P,'min_neighbors',0);
P = impose_default_value(P,'max_neighbors',1000);
P = impose_default_value(P,'indel_min_neighbors',1);
if P.indel_min_neighbors<1, error('indel_min_neighbors must be at least 1'); end
P = impose_default_value(P,'indel_max_neighbors',1000);
P = impose_default_value(P,'qual_min',0.1);
P = impose_default_value(P,'indel_qual_min',0.05);
P = impose_default_value(P,'min_territory_ratio',0);
P = impose_default_value(P,'indel_min_territory_ratio',10);
P = impose_default_value(P,'num_neighbor_patients',1);
P = impose_default_value(P,'impute_full_cov_when_promotes_significance',true);
P = impose_default_value(P,'restrict_to_one_mutation_per_patient',true);
P = impose_default_value(P,'radius_to_impute_coverage_around_mutations',10);
P = impose_default_value(P,'permutations_min_effect_size',1.01);
if P.permutations_min_effect_size<1, error('permutations_min_effect_size must be at least 1.00'); end
P = impose_default_value(P,'max_coverage_bins',10);
P = impose_default_value(P,'clustering_metric',204);
P = impose_default_value(P,'randseed',6789);
P = impose_default_value(P,'skip_permutations',false);
P = impose_default_value(P,'maxperm',1e5); %to speedup initial screens; need to reset to 1e6 later.
P = impose_default_value(P,'theta',1);
P = impose_default_value(P,'keyboard_before_begin',false);
P = impose_default_value(P,'penalty_per_strike',1.10);
if P.penalty_per_strike<1, error('P.penalty_per_strike is multiplicative, must be >=1'); end
P = impose_default_value(P,'extra_penalty_for_highest_noise_channels',0);
P = impose_default_value(P,'remove_duplicate_patients',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LOAD   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('LOADING DATA\n');

M=[];

% open FWB tracks
M.FWB = [];
M.FWB.conservation = org.broadinstitute.cga.tools.seq.FixedWidthBinary(conservation_fwb_file);
M.FWB.conservation.setNullVal(200);
M.FWB.context_and_effect = org.broadinstitute.cga.tools.seq.FixedWidthBinary(context_and_effect_fwb_file);
M.context_and_effect = load_struct(context_and_effect_categs_file);
M.context_and_effect.context65 = map_categories_to_65(context_and_effect_categs_file);
if ~strcmpi(basewise_coverage_fwb_file,'IMPUTE_FULL_COVERAGE')
  M.FWB.basewise_coverage = org.broadinstitute.cga.tools.seq.FixedWidthBinary(basewise_coverage_fwb_file);
  M.FWB.basewise_coverage.setNullVal(0);
end

% TARGET LIST
fprintf('Processing target list.\n');
M.targ = load_struct(target_list_file);
demand_fields(M.targ,{'gene','chr','start','end'});
M.targ.chr = convert_chr(M.targ.chr);
M.targ = make_numeric(M.targ,{'start','end'});
M.targ = reorder_struct_exclude(M.targ,isnan(M.targ.chr)|isnan(M.targ.start)|isnan(M.targ.end));
if any(M.targ.start>M.targ.end), error('target list has start>end'); end
M.targ = sort_struct(M.targ,{'gene','chr','start'});
M.targ.len = M.targ.end-M.targ.start+1;
if mean(M.targ.len)>5000, error('Looks like target list is whole transcripts.  Need individual exons!'); end

% MUTATIONS 
fprintf('Loading mutations...\n');
M.mut = load_struct(mutation_file);

% compute newbase if necessary
if isfield(M.mut,'newbase') && ~any(cellfun('isempty',M.mut.newbase))
  % already ok
elseif isfield(M.mut,'Reference_Allele') && isfield(M.mut,'Tumor_Seq_Allele1') && isfield(M.mut,'Tumor_Seq_Allele2')
  M.mut.newbase = M.mut.Tumor_Seq_Allele1;
  idx = find(strcmp(M.mut.Reference_Allele,M.mut.Tumor_Seq_Allele1));
  M.mut.newbase(idx) = M.mut.Tumor_Seq_Allele2(idx);
elseif isfield(M.mut,'ref_allele') && isfield(M.mut,'tum_allele1') && isfield(M.mut,'tum_allele2')
  M.mut.newbase = M.mut.tum_allele1;
  idx = find(strcmp(M.mut.ref_allele,M.mut.tum_allele1));
  M.mut.newbase(idx) = M.mut.tum_allele2(idx);
end

%%field processing
%if specified, explicitly use given fieldnames
%if isfield(P, 'fieldoverride'), %will clobber existing fields
%  demand_file(P.fieldoverride);
%  fo = load_struct_noheader(P.fieldoverride);  
%  M.mut = rename_fields(M.mut, fo.col1, fo.col2);
%end

%infer ambiguous fieldnames
%for now, this is just done for Variant_(classification|type) vs. (type|classification)
if P.infer_fields,
  D = load_struct(mutation_type_dictionary_file);
  cl = unique(D.classification); ty = unique(D.type);

  mfields = {'type', 'Variant_Classification', 'classification', 'Variant_Type'}; %possible fields
  ffields = {'classification', 'type'}; %final fields

  fnames = mfields(ismember(mfields, fieldnames(M.mut)));
  for i = 1:length(fnames),
    fn = fnames{i}; 
    [~, idx] = min([nnz(isnan(listmap(M.mut.(fn), cl))) nnz(isnan(listmap(M.mut.(fn), ty)))]);

    if isfield(M.mut, fn),
      M.mut = rename_field(M.mut, fn, ffields{idx});
    end
  end
end

% make sure required fields exist (renaming if necessary)
rf = {
    {'gene','Hugo_Symbol','Gene_name'},
    {'patient','Tumor_Sample_Barcode','Patient_name'},
    {'chr','Chromosome'},
    {'pos','Position','start','Start_position'},
    {'ref_allele','Reference_Allele'},
    {'newbase','Tumor_Allele','Tum_allele','Alt_allele','Alternate_allele','Tumor_Seq_Allele2'},
    {'type','Variant_Classification'},
    {'classification','Variant_Type'},
};

f = fieldnames(M.mut);
for i=1:length(rf)
  matches = find(ismember(lower(f),lower(rf{i})));
  if isempty(matches)
    fprintf('\nMutation file is missing a column named one of the following:\n');
    pr(rf{i});
    error('Mutation file missing %s column.',rf{i}{1});
  end
  if length(matches)>1
    fprintf('\nMutation file contains multiple columns for %s info:\n', rf{i}{1});
    pr(f(matches));
    match_to_first = find(strcmpi(lower(f),lower(rf{i}{1})));
    if ~isempty(match_to_first), matches = match_to_first; else matches = matches(1); end
    fprintf('Will use %s\n',f{matches});
  end
  M.mut = rename_field(M.mut,f{matches},rf{i}{1});
end

% hg18 liftover (if requested)
if strcmp(build, 'hg18'),
  disp('Performing hg18 liftover...');
  M.mut.pos18 = M.mut.pos;
  M.mut.pos = liftover(M.mut.chr, M.mut.pos, 'hg18', 'hg19');
  M.mut = reorder_struct(M.mut, ~isnan(M.mut.pos));
end

% remove duplicate mutations
[u ui uj] = unique_combos(M.mut.patient,M.mut.chr,M.mut.pos);
if length(ui)<slength(M.mut)
  fprintf('Keeping %d/%d unique mutations.\n',length(ui),slength(M.mut));
  M.mut = reorder_struct(M.mut,ui);
end

% remove duplicate patients
fprintf('Scanning for duplicate patients...\n');
X = new_find_duplicate_samples(M.mut);
if ~isempty(X.drop)
  fprintf('Removing the following %d duplicate patients:\n',length(X.drop));
  disp(X.drop);
  M.mut = reorder_struct_exclude(M.mut,ismember(M.mut.patient,X.drop));
end

% apply blacklist (if specified)
if ~strcmp('',mutation_blacklist_file) && ~strcmp('none',mutation_blacklist_file)
  fprintf('Applying blacklist: %s\n',mutation_blacklist_file);
  M.mut = apply_mutation_blacklist(M.mut,mutation_blacklist_file,P);
end

% COVERAGE 

fprintf('Loading coverage models ... ');
load(coverage_models_mat_file,'C');
demand_fields(C,{'type','ntype','gene','cat','ncat','gene_effect_terr','gene_effect_cov'});
if C.ncat~=1885, error('wrong C.ncat'); end
if any(isinf(C.gene_effect_terr(:))) || any(isinf(C.gene_effect_cov(:))), error('infs in terr/cov'); end
if any(isnan(C.gene_effect_terr(:))) || any(isnan(C.gene_effect_cov(:))), error('nans in terr/cov'); end
M.cov = C; clear C;
fprintf('%d scheme(s).\n', M.cov.ntype);

% enforce target list (if requested)
if P.enforce_target_list
  fprintf('Enforcing target list.\n');
  PP=P;
  PP.enforce_target_list_footprint_table = M.cov.gene;
  tmp=[]; tmp.map = map_mutations_to_targets(M.mut,M.targ,PP);
  idx = find(isnan(tmp.map));
  if ~isempty(idx)
    fprintf('Removing %d/%d mutations that fall outside target gene intervals.\n',length(idx),slength(M.mut));
    M.mut = reorder_struct_exclude(M.mut,idx);
    tmp.map(idx)=[];
  end
  tmp.new_gene = M.targ.gene(tmp.map);
  if isfield(M.mut,'gene')
    tmp.old_gene = M.mut.gene;
    idx = find(~strcmpi(tmp.old_gene,tmp.new_gene));
    if ~isempty(idx)
      fprintf('Reassigning the following %d gene identities:', length(idx));
      tmp.reassign = stringsplice([tmp.old_gene(idx) tmp.new_gene(idx)],1,' -> ');
      count(tmp.reassign,1);
    end
  end
  M.mut.gene = tmp.new_gene;
end
tmp = [];

% remove IGR mutations
idx = find(strcmp('',M.mut.gene)|strcmpi('unknown',M.mut.gene));
if ~isempty(idx)
  fprintf('Omitting %d/%d mutations that are intergenic or unknown gene.\n',length(idx),slength(M.mut));
  M.mut = reorder_struct_exclude(M.mut,idx);
end

% EFFECT 
M.effect = [];
M.effect.name = {'ncd';'syn';'mis';'non';'spl';'indel_ncd';'indel_cod';'indel_spl'};
M.effect.permutations_effect_idx = [1;2;3;4;4;1;5;5];
M.effect.include_in_permutations = (M.effect.permutations_effect_idx>=3);
M.effect.is_indel = grepmi('indel',M.effect.name);

fprintf('Looking up "effect" in mutation_type_dictionary_file\n');
D = load_struct(mutation_type_dictionary_file);
demand_fields(D,'effect');
f = setdiff(fieldnames(D),{'effect'});
demand_fields(M.mut,f);
D.effect_idx = listmap(D.effect,M.effect.name);
di = multimap(M.mut,D,f);
M.mut.effect_idx = nansub(D.effect_idx,di);
idx = find(isnan(M.mut.effect_idx));
if ~isempty(idx)
  fprintf('Omitting %d/%d mutations of unknown "effect"\n',length(idx),slength(M.mut));
  M.mut = reorder_struct_exclude(M.mut,idx);
end

% other mutation conversions
fprintf('Converting mutation data...\n');
M.mut.chr = convert_chr(M.mut.chr);
M.mut = make_numeric(M.mut,'pos');
M.mut.newbase_idx = listmap(regexprep(M.mut.newbase,'^(.).*$','$1'), {'A', 'C', 'G', 'T'});
M.mut.newbase_idx(M.effect.is_indel(M.mut.effect_idx)) = 5;  % indels get newbase_idx=5
idx = find(isnan(M.mut.newbase_idx));
if ~isempty(idx)
  fprintf('Omitting %d/%d mutations of invalid "newbase"\n',length(idx),slength(M.mut));
  M.mut = reorder_struct_exclude(M.mut,idx);
end
M.mut = reorder_struct(M.mut, ~isnan(M.mut.chr));

% PATIENTS 

[M.pat.name tmp M.mut.pat_idx] = unique(M.mut.patient);
M.np = slength(M.pat);
fprintf('%d patients\n',M.np);
if M.np<2, fprintf('WARNING:  MutSig is not applicable to single patients.\n'); end

% impute callschemes
if M.cov.ntype == 1,
  fprintf('Using single coverage model "%s" for all patients.\n', M.cov.type{1})
  M.pat.callscheme_name = repmat(M.cov.type(1), slength(M.pat), 1);
  M.pat.callscheme = zeros(slength(M.pat), 1);

  M.mut = add_and_convert_simple_fieldnames(M.mut);
  M.mut = add_helper_is_fields(M.mut); 

  pvec = 1:slength(M.pat); 
  M.pat.nmut = as_column(histc(M.mut.pat_idx,1:slength(M.pat)));
  M.pat.n_coding = as_column(histc(M.mut.pat_idx(M.mut.is_coding), pvec));
  M.pat.n_flank = as_column(histc(M.mut.pat_idx(M.mut.is_flank), pvec));
  M.pat.n_coding_nonsilent = as_column(histc(M.mut.pat_idx(M.mut.is_coding & ~M.mut.is_silent), pvec));
  M.pat.fracflank = M.pat.n_flank./M.pat.nmut;
else
  M = impute_callschemes(M);
end

if ~all(M.pat.callscheme>=0 & M.pat.callscheme<=3), error('Problem imputing callschemes'); end
M.pat.cov_idx = mapacross(M.pat.callscheme,[0 1 2 3],[1 2 3 5]);
M.pat = rmfield(M.pat,'callscheme');

% remove genes with very low coverage of coding regions
% (based on the coverage models that were identified as being relevant)
cod = grepv('noncoding',M.cov.cat.name,1);
gene_terr = sum(M.cov.gene_effect_terr(:,:,cod),3);
gene_model_cov = sum(M.cov.gene_effect_cov(:,:,cod),3);
num_each_model = histc(M.pat.cov_idx,(1:M.cov.ntype)');
gene_tot_terr = gene_terr * M.np;
gene_tot_cov = gene_model_cov * num_each_model;
gene_frac_cov = gene_tot_cov ./ gene_tot_terr;

idx = find(gene_frac_cov<P.gene_min_frac_coverage_required);
if ~isempty(idx)
  oldng = M.cov.ng;
  fprintf('Omitting %d/%d genes because they have extremely low coverage.\n',length(idx),oldng);
  M.cov.gene = reorder_struct_exclude(M.cov.gene,idx);
  M.cov.ng = M.cov.ng - length(idx);
  M.cov.gene_effect_cov(idx,:,:)=[];
  M.cov.gene_effect_terr(idx,:,:)=[];
end

% GENES 

M.gene = M.cov.gene; M.cov = rmfield(M.cov,'gene');
M.ng = slength(M.gene); M.cov = rmfield_if_exist(M.cov,'ng');
M.targ.gene_idx = listmap(M.targ.gene,M.gene.name);
M.mut.gene_idx = listmap(M.mut.gene,M.gene.name);
idx = find(isnan(M.mut.gene_idx));
if ~isempty(idx)
  fprintf('Omitting %d/%d mutations outside gene set.\n',length(idx),slength(M.mut));
  M.mut = reorder_struct_exclude(M.mut,idx);
end

% CATEGORIES 

% look up context65 if necessary
if isfield(M.mut,'context65')
  M.mut = make_numeric(M.mut,'context65');
else
  fprintf('\nLooking up values in context_and_effect track\n');
  M.mut.context_and_effect = double(M.FWB.context_and_effect.get(M.mut.chr,M.mut.pos));
  M.mut.context_and_effect(M.mut.context_and_effect==-1)=nan;
  M.mut.context65 = nansub(M.context_and_effect.context65,M.mut.context_and_effect);
end

% collapse coverage models, first step: collapse 1885 to 192x5degree
M.cov.gene_effect_terr = collapse_1885_to_192x5(M.cov.gene_effect_terr,3);
M.cov.gene_effect_cov = collapse_1885_to_192x5(M.cov.gene_effect_cov,3);
M.cov.dims_of_gene_effect_cov = 'gene, covmodel, 192categories, effect (ncd/syn/mis/non/spl)';

fprintf('Category discovery....\n');

M.effect.include_in_category_discovery = grepm('^(syn|mis|non|spl)$',M.effect.name);
midx = find(nansub(M.effect.include_in_category_discovery,M.mut.effect_idx)==1);
n = hist2d(M.mut.context65(midx),M.mut.newbase_idx(midx),1:65,1:4);
N = squeeze(sum(sum(M.cov.gene_effect_cov(:,:,:,[2:end]),4),1))';
if size(N, 1) ~= 192, N = N'; end %in the case that we only have one callscheme
N = 3*N(1:3:end,:);
N = round(N * histc(M.pat.cov_idx,1:M.cov.ntype));
N(65)=0;
Nn = collapse_Nn_65_to_32([N n]);
PP=[]; PP.max_k = P.number_of_categories_to_discover;
PP.mutcategs_report_filename = [output_dir '/mutcateg_discovery.txt'];
Ks = find_mut_categs(Nn,PP);
M.categ = Ks{P.number_of_categories_to_discover};
M.ncat = slength(M.categ);

% collapse coverage models, second step: collapse 192x5degree to ncatx5degree
M.cov.gene_effect_terr = collapse_192_to_categ_set(M.cov.gene_effect_terr,M.categ,3);
M.cov.gene_effect_cov = collapse_192_to_categ_set(M.cov.gene_effect_cov,M.categ,3);
M.cov.dims_of_gene_effect_cov = 'gene, covmodel, ncat, effect (ncd/syn/mis/non/spl)';
M.cov = move_field_to_before(M.cov,'dims_of_gene_effect_cov','gene_effect_cov');
M.cov = rmfield(M.cov,{'cat','ncat'}); % remove old category information

% analyze collapsed categories vs. full-spectrum context_and_effect categories
M.context_and_effect.Q = collapse_context_and_effect_categories(M.context_and_effect,M.categ);

% assign mutations to categories
fprintf('Assigning mutation categories...\n');
M.mut = rmfield_if_exist(M.mut,{'categ','categ_ignoring_null_categ'});
c = assign_65x4_to_categ_set(M.categ);
c = bsxfun(@times,c,1:slength(M.categ));
c = squeeze(max(c,[],2));
M.mut.categ_idx = zeros(slength(M.mut),1);
for newbase_idx=1:4
  idx = find(M.mut.newbase_idx==newbase_idx);
  M.mut.categ_idx(idx) = nansub(c(:,newbase_idx),M.mut.context65(idx));
end
idx = find(M.mut.categ_idx==0 & ~M.effect.is_indel(M.mut.effect_idx));
if ~isempty(idx)
  fprintf('Omitting %d/%d mutations of unassigned categ.\n',length(idx),slength(M.mut));
  M.mut = reorder_struct_exclude(M.mut,idx);
end

% COVARIATES 

fprintf('Loading covariate files.\n');
V = load_struct_specify_string_cols(covariates_file,1);  % gene is string
f = fieldnames(V); if ~strcmp(f{1},'gene'), error('first column of covariates file must be "gene"'); end
M.cvname = f(2:end); M.nv = length(M.cvname);
gidx = listmap(M.gene.name,V.gene);
M.V = nan(M.ng,M.nv);
for vi=1:M.nv, M.V(:,vi) = nansub(V.(M.cvname{vi}),gidx); end

% convert covariate raw values to Z-scores
M.Z = nan(M.ng,M.nv);
for vi=1:M.nv
  missing = isnan(M.V(:,vi)) | isinf(M.V(:,vi));
  mn = mean(M.V(~missing,vi));
  sd = std(M.V(~missing,vi),0);  % second parameter=0 means confirm default behavior of normalize by (N-1) not (N)
  M.Z(~missing,vi) = (M.V(~missing,vi)-mn)./sd;
end

% compute number of "strikes"
if all(isfield(M.gene,{'log_exprmax','rt','hiC','paz'}))
  % legacy method
  M.gene.nstrikes = (M.gene.log_exprmax<4.5)+(M.gene.log_exprmax<4.0)+(M.gene.rt>600)+(M.gene.rt>800)+...
      (M.gene.hiC<-0.02)+(M.gene.hiC<-0.01)+(M.gene.paz>0.2)+(M.gene.paz>0.3);
else
  % approximate general method for future release
  % (needs to be replaced with a bootstrapping method)
  M.gene.nstrikes = sum(M.Z>1.5,2) + sum(M.Z>2,2);
end

% compute per-gene minimum effect size
%M.gene.min_effect_size = 1.01*1.10.^(M.gene.nstrikes);
M.gene.min_effect_size = 1.01*P.penalty_per_strike.^(M.gene.nstrikes);

% compute geneidx lookup table for speedup
M.mut = sort_struct(M.mut,'gene_idx');
[u ui uj] = unique(M.mut.gene_idx,'first');
h = histc(uj,1:length(u));
M.geneidx = cell(M.ng,1);
for i=1:length(u), M.geneidx{u(i)} = as_column(ui(i):ui(i)+h(i)-1);end

% increase territory to the max coverage represented
cov_idx = unique(M.pat.cov_idx);
M.cov.gene_effect_terr = ceil(max([M.cov.gene_effect_terr M.cov.gene_effect_cov(:,cov_idx,:,:)],[],2));

% TOTAL COUNTS

M.gene.codelen = round(sum(sum(M.cov.gene_effect_terr(:,1,:,2:5),4),3)/3);

M.gene.Nncd=0;
M.gene.Nsyn=0;
M.gene.Nmis=0;
M.gene.Nnon=0;
M.gene.Nspl=0;
for cov_idx=1:M.cov.ntype
  pidx = find(M.pat.cov_idx==cov_idx);
  M.gene.Nncd = M.gene.Nncd + sum(M.cov.gene_effect_cov(:,cov_idx,:,1),3)*length(pidx);
  M.gene.Nsyn = M.gene.Nsyn + sum(M.cov.gene_effect_cov(:,cov_idx,:,2),3)*length(pidx);
  M.gene.Nmis = M.gene.Nmis + sum(M.cov.gene_effect_cov(:,cov_idx,:,3),3)*length(pidx);
  M.gene.Nnon = M.gene.Nnon + sum(M.cov.gene_effect_cov(:,cov_idx,:,4),3)*length(pidx);
  M.gene.Nspl = M.gene.Nspl + sum(M.cov.gene_effect_cov(:,cov_idx,:,5),3)*length(pidx);
end
M.gene.Nind = M.gene.Nncd + M.gene.Nsyn + M.gene.Nmis + M.gene.Nnon + M.gene.Nspl;

M.gene.Nncd = round(M.gene.Nncd);
M.gene.Nsyn = round(M.gene.Nsyn);
M.gene.Nmis = round(M.gene.Nmis);
M.gene.Nnon = round(M.gene.Nnon);
M.gene.Nspl = round(M.gene.Nspl);
M.gene.Nind = round(M.gene.Nind);

M.gene.nncd = as_column(histc(M.mut.gene_idx(M.mut.effect_idx==find(strcmp('ncd',M.effect.name))),1:M.ng));
M.gene.nsyn = as_column(histc(M.mut.gene_idx(M.mut.effect_idx==find(strcmp('syn',M.effect.name))),1:M.ng));
M.gene.nmis = as_column(histc(M.mut.gene_idx(M.mut.effect_idx==find(strcmp('mis',M.effect.name))),1:M.ng));
M.gene.nnon = as_column(histc(M.mut.gene_idx(M.mut.effect_idx==find(strcmp('non',M.effect.name))),1:M.ng));
M.gene.nspl = as_column(histc(M.mut.gene_idx(M.mut.effect_idx==find(strcmp('spl',M.effect.name))),1:M.ng));
M.gene.nind = as_column(histc(M.mut.gene_idx(ismember(M.mut.effect_idx,find(grepmi('indel_cod|indel_spl',M.effect.name)))),1:M.ng));

% compute gene "f" (gene-specific factor) from its own synonymous mutations
M.gene.f = (M.gene.nsyn ./ M.gene.Nsyn) / (sum(M.gene.nsyn)/sum(M.gene.Nsyn));

% TOTAL RATES 

globalrate_ncd = sum(M.gene.nncd) / sum(M.gene.Nncd);
globalrate_syn = sum(M.gene.nsyn) / sum(M.gene.Nsyn);
globalrate_mis = sum(M.gene.nmis) / sum(M.gene.Nmis);
globalrate_non = sum(M.gene.nnon) / sum(M.gene.Nnon);
globalrate_spl = sum(M.gene.nspl) / sum(M.gene.Nspl);
globalrate_ind = sum(M.gene.nind) / sum(M.gene.Nind);

fprintf('Global rates (/Mb):  ncd %.2f   syn %.2f   mis %.2f  non %.2f  spl %.2f  ind %.2f\n',...
   globalrate_ncd*1e6,globalrate_syn*1e6,globalrate_mis*1e6,globalrate_non*1e6,globalrate_spl*1e6,globalrate_ind*1e6);
if globalrate_mis==0, error('zero missense rate'); end

% PATIENT RATES 

M.pat.N_tot = nan(M.np,1);       % total noncoding+coding territory (for total patient rate calculation)
M.pat.N_c = nan(M.np,M.ncat);    % total coding territory broken down by context category (for per-category calculation)
for c=1:M.cov.ntype
  pidx = find(M.pat.cov_idx==c);
  M.pat.N_tot(pidx) = fullsum(M.cov.gene_effect_cov(:,c,:,:));
  M.pat.N_c(pidx,:) = repmat(as_row(squeeze(sum(sum(M.cov.gene_effect_cov(:,c,:,2:5),1),4))),length(pidx),1);
end
M.pat.N_tot = round(M.pat.N_tot/3);
M.pat.N_c = round(M.pat.N_c/3);
M.pat.nsil_tot = histc(M.mut.pat_idx(M.mut.effect_idx==2),1:M.np);
M.pat.nnon_tot = histc(M.mut.pat_idx(ismember(M.mut.effect_idx,[3:5])),1:M.np);
M.pat.n_tot = histc(M.mut.pat_idx,1:M.np);
M.pat.rate_sil = M.pat.nsil_tot./M.pat.N_tot;
M.pat.rate_non = M.pat.nnon_tot./M.pat.N_tot;
M.pat.rate_tot = M.pat.n_tot./M.pat.N_tot;
M.pat.log_rate_tot = max(-9,log10(M.pat.rate_tot));

% also compute table of (patient,category)-specific rates
midx = find(ismember(M.mut.effect_idx,find(ismember(M.effect.name,{'syn','mis','non','spl'}))));
M.pat.n_c = hist2d_fast(M.mut.pat_idx(midx),M.mut.categ_idx(midx),1,M.np,1,M.ncat);
midx = find(ismember(M.mut.effect_idx,find(ismember(M.effect.name,{'indel_cod','indel_spl'}))));
M.pat.n_ind = as_column(histc(M.mut.pat_idx(midx),1:M.np));
M.pat.N_ind = sum(M.pat.N_c,2);
M.pat.rate_c = M.pat.n_c ./ M.pat.N_c;
M.pat.rate_ind = M.pat.n_ind ./ M.pat.N_ind;

% NEIGHBOR PATIENTS
if P.num_neighbor_patients>0
  % compute matrix of neighbor-patients
  if P.num_neighbor_patients>=20
    patneighb = zeros(M.np,M.np);
  else
    patneighb = sparse(M.np,M.np);
  end
  for p=1:M.np
    patdist = abs(M.pat.log_rate_tot-M.pat.log_rate_tot(p));
    [tmp ord] = sort(patdist);
    npidx = ord(1:min(length(ord),P.num_neighbor_patients+1));   % (include the patient itself)
    if isfield(M.pat,'cov_idx')
      % restrict neighbor patients to within the same callscheme (avoids an N=0 bug)
      npidx = intersect(npidx,find(M.pat.cov_idx==M.pat.cov_idx(p)));
    end
    patneighb(p,npidx)=1;
  end
end

% will discover two sets of bagels:
%      main_bagel
%      indel_bagel

n_bageltypes = 2;
min_neighbors = [P.min_neighbors P.indel_min_neighbors];
max_neighbors = [P.max_neighbors P.indel_max_neighbors];
min_territory_ratio = [P.min_territory_ratio P.indel_min_territory_ratio];
qual_min = [P.qual_min P.indel_qual_min];

max_min_neighbors = max(min_neighbors);
max_max_neighbors = max(max_neighbors);
max_min_territory_ratio = max(min_territory_ratio);
min_qual_min = min(qual_min);

M.gene.bagel = nan(M.ng,max_max_neighbors);
M.gene.nnei = zeros(M.ng,1);
M.gene.nnei2 = M.gene.nnei;

% PERMUTATIONS preparation

% initialize random number generator
fprintf('Permutations: initializing with randseed=%d\n',P.randseed);
rand('twister',P.randseed);

% field for "position along track" that will be filled out during procedure
M.mut.trackpos = nan(slength(M.mut),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% RUN 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nProcessing each gene...\n');

z = nan(M.ng,1);
M.gene.npat = z;
M.gene.nsite = z;
M.gene.pCVmin = z;
M.gene.pCVmid = z;
M.gene.pCVmax = z;
M.gene.nmut = z;
M.gene.nperm = z;
M.gene.pCL = z;
M.gene.pFN = z;
M.gene.pCLFN = z;

tt=tic();
qqhist = zeros(1,6);

fprintf('now       eta                # gene        nmut   nperm  CV          CL          FN           [ qq plot ]\n');

if P.keyboard_before_begin, keyboard; end

for g=1:M.ng
 gene_status_reported=false;
 pCVmax = 1; pclust=1; pcons=1; nperm=0; nm=0;
 try

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % BAGEL FINDING
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % calculate distances from this gene
  df2 = bsxfun(@minus,M.Z,M.Z(g,:)).^2;
  dist2 = nansum(df2,2)./sum(~isnan(df2),2);
  [tmp,ord] = sort(dist2); ord = [g;ord(ord~=g)];

  per_type_nnei = zeros(n_bageltypes,1);
  per_type_bagel = repmat({nan(1,max_max_neighbors)},n_bageltypes,1);
  per_type_bagel_done = false(n_bageltypes,1);

  % expand bagel outward until quality falls below qual_min
  nfit=0; Nfit=0;
  for ni=0:max_max_neighbors, gidx = ord(ni+1);

    ngene = M.gene.nncd(gidx) + M.gene.nsyn(gidx);
    Ngene = M.gene.Nncd(gidx) + M.gene.Nsyn(gidx);

    if ni==0, ngene0=ngene; Ngene0=Ngene; end
    nfit=nfit+ngene; Nfit=Nfit+Ngene;

    %hack to fix bug with coverage model:
    if Ngene < ngene || Ngene0 < ngene0, continue; end

    % compare the gene being added to the central gene
    qual = 2*hyge2cdf(ngene,Ngene,ngene0,Ngene0);  % (two-sided)
    if qual>1, qual = 2-qual; end

    % territory ratio
    territory_ratio = Nfit/Ngene0;

    % stopping criterion: stop if this gene would drop quality below qual_min
    if ((ni-1)>=max_min_neighbors && territory_ratio>=max_min_territory_ratio) && qual<min_qual_min, break; end

    % update gene's bagel
    if ni>0
      M.gene.bagel(g,ni)=gidx;   % NOTE: not actually used anywhere
      % which categories include this bagel?
      for c=1:n_bageltypes
        if ~per_type_bagel_done(c) && (ni<=min_neighbors(c) || territory_ratio<min_territory_ratio(c) || (ni<=max_neighbors(c) && qual>=qual_min(c)))
          per_type_nnei(c) = per_type_nnei(c) + 1;
          per_type_bagel{c}(ni) = gidx;
        else
          per_type_bagel_done(c) = true;
        end
      end
    end
  end % next neighborhood size

  %remove entries 
  for c = 1:n_bageltypes,
    per_type_bagel{c}((per_type_nnei(c) + 1):end) = [];
    per_type_bagel{c}(isnan(per_type_bagel{c})) = [];
  end

  M.gene.nnei(g) = per_type_nnei(1);   % main bagel
  M.gene.nnei2(g) = per_type_nnei(2);  % indel bagel

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PROJECTION METHOD
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % STEP 0
  % retrieve counts for this gene + bagel

  % gene's mutation counts
  gene_midx = M.geneidx{g};
  midx1 = gene_midx(M.mut.effect_idx(gene_midx)==find(strcmp('mis',M.effect.name)));
  gene_nmis = hist2d_fast(M.mut.pat_idx(midx1),M.mut.categ_idx(midx1),1,M.np,1,M.ncat);
  midx2 = gene_midx(M.mut.effect_idx(gene_midx)==find(strcmp('non',M.effect.name)));
  gene_nnon = hist2d_fast(M.mut.pat_idx(midx2),M.mut.categ_idx(midx2),1,M.np,1,M.ncat);
  midx3 = gene_midx(M.mut.effect_idx(gene_midx)==find(strcmp('spl',M.effect.name)));
  gene_nspl = hist2d_fast(M.mut.pat_idx(midx3),M.mut.categ_idx(midx3),1,M.np,1,M.ncat);
  midx4 = gene_midx(ismember(M.mut.effect_idx(gene_midx),find(grepm('indel_(cod|spl)',M.effect.name))));
  if ~isempty(midx4), gene_nind = as_column(histc(M.mut.pat_idx(midx4),1:M.np)); else gene_nind = zeros(M.np,1); end

  % find # unique patients and sites (only for display in table; not used directly in calculation)
  midx = [midx1;midx2;midx3;midx4];
  if isempty(midx)
    M.gene.nsite(g) = 0;
    M.gene.npat(g) = 0;
  else
    M.gene.nsite(g) = length(unique_combos(M.mut.chr(midx),M.mut.pos(midx)));
    M.gene.npat(g) = length(unique(M.mut.pat_idx(midx)));
  end

  % gene's coverage (or territory)
  if P.impute_full_cov_when_promotes_significance
    % use territory
    tmp = squeeze(M.cov.gene_effect_terr(g,1,:,2:5));
    gene_Nsyn = repmat(tmp(:,1)',M.np,1);
    gene_Nmis = repmat(tmp(:,2)',M.np,1);
    gene_Nnon = repmat(tmp(:,3)',M.np,1);
    gene_Nspl = repmat(tmp(:,4)',M.np,1);
  else  % use coverage
    gene_Nsyn = zeros(M.np,M.ncat);
    gene_Nmis = zeros(M.np,M.ncat);
    gene_Nnon = zeros(M.np,M.ncat);
    gene_Nspl = zeros(M.np,M.ncat);
    for c=1:M.ncat
      gene_Nsyn(:,c) = squeeze(M.cov.gene_effect_cov(g,M.pat.cov_idx,c,2));
      gene_Nmis(:,c) = squeeze(M.cov.gene_effect_cov(g,M.pat.cov_idx,c,3));
      gene_Nnon(:,c) = squeeze(M.cov.gene_effect_cov(g,M.pat.cov_idx,c,4));
      gene_Nspl(:,c) = squeeze(M.cov.gene_effect_cov(g,M.pat.cov_idx,c,5));
    end
  end
  gene_Nind = sum(gene_Nsyn + gene_Nmis + gene_Nnon + gene_Nspl,2);

  % main sphere (=gene+bagel): noncoding and synonymous
  main_sphere_nncd = zeros(M.np,M.ncat);
  main_sphere_Nncd = zeros(M.np,M.ncat);
  main_sphere_nsyn = zeros(M.np,M.ncat);
  main_sphere_Nsyn = zeros(M.np,M.ncat);
  bageltype = 1;
  main_sphere = [g per_type_bagel{bageltype}];
  main_sphere_midx = cat(1,M.geneidx{main_sphere});
  for c=1:M.ncat
    midx = main_sphere_midx(M.mut.effect_idx(main_sphere_midx)==find(strcmp('ncd',M.effect.name)) & M.mut.categ_idx(main_sphere_midx)==c);
    if ~isempty(midx), main_sphere_nncd(:,c) = as_column(histc(M.mut.pat_idx(midx),1:M.np)); end
    midx = main_sphere_midx(M.mut.effect_idx(main_sphere_midx)==find(strcmp('syn',M.effect.name)) & M.mut.categ_idx(main_sphere_midx)==c);
    if ~isempty(midx), main_sphere_nsyn(:,c) = as_column(histc(M.mut.pat_idx(midx),1:M.np)); end
    for cov_idx=1:M.cov.ntype
      pidx = (M.pat.cov_idx==cov_idx);
      main_sphere_Nncd(pidx,c) = sum(M.cov.gene_effect_cov(main_sphere,cov_idx,c,1),1);
      main_sphere_Nsyn(pidx,c) = sum(M.cov.gene_effect_cov(main_sphere,cov_idx,c,2),1);
    end
  end
  % expand to neighboring patients (if specified)
  if P.num_neighbor_patients>0
    main_sphere_nncd = patneighb*main_sphere_nncd;
    main_sphere_Nncd = patneighb*main_sphere_Nncd;
    main_sphere_nsyn = patneighb*main_sphere_nsyn;
    main_sphere_Nsyn = patneighb*main_sphere_Nsyn;
  end

  % indel bagel
  indel_bagel_nind = zeros(M.np,1);
  indel_bagel_Nind = zeros(M.np,1);
  bageltype = 2;
  indel_bagel = per_type_bagel{bageltype};
  if isempty(indel_bagel),
    indel_bagel_Nind = 0;
  else
    indel_bagel_midx = cat(1,M.geneidx{indel_bagel});
    midx = indel_bagel_midx(ismember(M.mut.effect_idx(indel_bagel_midx),find(grepm('indel_(cod|spl)',M.effect.name))));
    if ~isempty(midx), indel_bagel_nind = as_column(histc(M.mut.pat_idx(midx),1:M.np)); end
    for cov_idx=1:M.cov.ntype
      pidx = (M.pat.cov_idx==cov_idx);
      indel_bagel_Nind(pidx) = sum(sum(sum(M.cov.gene_effect_cov(indel_bagel,cov_idx,:,:),3),4),1);
    end
  end

  % STEP 1
  % compute probability of seeing zero mutations in each degree

  main_ntot = gene_nmis + gene_nnon + gene_nspl + main_sphere_nncd + main_sphere_nsyn;
  main_Ntot = gene_Nmis + gene_Nnon + gene_Nspl + main_sphere_Nncd + main_sphere_Nsyn;
  indel_ntot = gene_nind + indel_bagel_nind;
  indel_Ntot = gene_Nind + indel_bagel_Nind;

  % extra penalty for highest-noise channels
  if P.extra_penalty_for_highest_noise_channels>0
    a = 1.1.^P.extra_penalty_for_highest_noise_channels; % a=0: no effect,   a=1, 1.1x for every 10xBMR    a=4, 1.4x for every 10xBMR
    w = @(x) a.^max(-4,min(4,log10(x/median(x(:)))));
    main_Ntot = main_Ntot ./ w(M.pat.rate_c);
    indel_Ntot = indel_Ntot ./ w(M.pat.rate_ind);
  end

  nm = sum((any((gene_nmis+gene_nnon+gene_nspl)>0,2)+gene_nind)>0);  % (only used for display in status report line)
  % TO DO:  to speed up runs on small datasets, add early escape for genes with no mutations

              %  1   2   3   4
  ndeg = 4;   % mis non spl ind
  P0 = nan(M.np,ndeg);
  has_mutation = false(M.np,ndeg);
  for d=1:ndeg
    if d<4
      n_total = main_ntot; N_total = main_Ntot;
      if d==1
        n_signal = gene_nmis; N_signal = gene_Nmis;
      elseif d==2
        n_signal = gene_nnon; N_signal = gene_Nnon;
      else % d==3
        n_signal = gene_nspl; N_signal = gene_Nspl;
      end
    else % d==4
      n_total = indel_ntot; N_total = indel_Ntot;
      n_signal = gene_nind; N_signal = gene_Nind;
    end
    % make sure we never have a more than 1000x difference between Nsignal and Ntotal: is not realistic
    N_total(N_total>1000*N_signal) = 1000*N_signal(N_total>1000*N_signal);     % to catch an N=0 bug
    % make sure values are rounded
    N_total = round(N_total); n_total = round(n_total); N_signal = round(N_signal);
    % make sure we never have N<2*n
    N_signal(2*n_signal>N_signal)=2*n_signal(2*n_signal>N_signal);
    N_total(2*n_total>N_total)=2*n_total(2*n_total>N_total);
    % find probabilities
    p = my_hygepdf(0,N_total,n_total,N_signal);  % (2D vectorization was being screwed up by all the error checking in hygepdf.m)
    P0(:,d) = prod(p,2);
    has_mutation(:,d) = sum(n_signal,2)>0;
  end
  P0(isnan(P0))=1;
  P0(P0>1)=1;
  P0(P0<0)=0;
  P1 = 1-P0;

  % for each sample, prioritize mutation categories according to how likely
  % it would be for this gene x sample to have a mutation there by chance.
  %
  % determine each patient's priority order of categories (according to P1)
  %  left column of "priority" =  lowest priority =  most likely category of mutation
  % right column of "priority" = highest priority = least likely category of mutation
  [tmp priority] = sort(P1,2,'descend');
  % sort the P arrays to put the columns in least->most priority order
  shft = (priority - repmat(1:ndeg,M.np,1));
  map = reshape(1:(M.np*ndeg),M.np,ndeg);
  newmap = map + shft*M.np;
  P0 = P0(newmap);
  P1 = P1(newmap);
  has_mutation = has_mutation(newmap);

  % STEP 2
  % for each sample, compute probability that it would have been of each degree.
  % where d=0 (no mut) ..... ndeg (most extreme mut)
  Pdeg = nan(M.np,ndeg+1);
  for d=0:ndeg
    % has to be clear in all categories > d
    Pdeg(:,d+1) = prod(P0(:,d+1:end),2);
    % and nonclear in category d (if d>0)
    if d>0, Pdeg(:,d+1) = Pdeg(:,d+1) .* P1(:,d); end
  end

  %% STEP 2a: calculate score for a sample being of each possible degree
  %% (uses new style, where score = -log10 probability of having a mutation in that category
  %% (zero score for no mutations)
  Sdeg = [zeros(M.np,1) -log10(P1)];

  % score multiplier (to increase resolution of scores)
  Sdeg = Sdeg * 10;

  % round all scores UP to nearest integer
  Sdeg = ceil(Sdeg);

  % score cap
  Sdeg = min(Sdeg,1000);   % because inf will crash the C version!

  % STEP 3
  % determine actual degree of each sample, and score_obs for the gene

  degree = zeros(M.np,1);
  score_obs = 0;
  for p = 1:M.np
    for d = ndeg:-1:1
      if has_mutation(p,d), degree(p) = d; break; end
    end
    score_obs = score_obs + Sdeg(p,degree(p)+1);
  end

  % impose minimum effect size by decreasing score_obs
  score_obs = score_obs / M.gene.min_effect_size(g);

  % STEP 4
  % compute P value for gene by convolutions

  numbins = ceil(score_obs+max(5,0.2*score_obs));
  H = zeros(numbins,1);
  newH = zeros(numbins,ndeg+1);

  [pCVmax pCVmin] = projection_1d_convolutions_fast(Sdeg,Pdeg,score_obs,numbins,H,newH);
  % pmax is the p-value generated by the standard procedure we've always used.
  % pmin is the (better) p-value that goes one discrete step further
  % to solve the problem of discrete statistics, we want to take a randomly chosen position between pmin and pmax.
  M.gene.pCVmax(g) = pCVmax;                       % for sorting genelist
  M.gene.pCVmid(g) = pCVmin+rand*(pCVmax-pCVmin);  % for q-q plot
  M.gene.pCVmin(g) = pCVmin;                       % for future reference

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PREPARATION FOR PERMUTATIONS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  while(true) % container loop, will only be executed once, but we can "break" out at any point to abort going into the permutations

    if P.skip_permutations || (P.maxperm==0), break; end

    % find the target regions for this gene
    tidx = find(M.targ.gene_idx==g);
    genelength = sum(M.targ.len(tidx));
    if genelength==0, break; end        % no targets found
  
    % find the mutations for this gene and map to targets
    midx = find(M.mut.gene_idx==g & M.effect.include_in_permutations(M.mut.effect_idx));
    exonstart=1;
    for ti=1:length(tidx),t=tidx(ti);
      z = find(M.mut.chr(midx)==M.targ.chr(t) & M.mut.pos(midx)>=M.targ.start(t) & M.mut.pos(midx)<=M.targ.end(t));
      M.mut.trackpos(midx(z)) = exonstart + M.mut.pos(midx(z))-M.targ.start(t);
      exonstart=exonstart+M.targ.len(t);
    end
    midx(isnan(M.mut.trackpos(midx)))=[];   % remove mutations that didn't map to a target
  
    if P.restrict_to_one_mutation_per_patient   % (chosen randomly)
      midx = midx(randperm(length(midx))); [u ui uj] = unique(M.mut.pat_idx(midx)); midx = midx(ui);
    end
    nm = length(midx);
    M.gene.nmut(g) = nm;
    if nm<2, break; end  % we only do permutations on genes with >=2 mutations
  
    % read conservation, coverage, and context_and_effect for these regions
    conservation_track = double(M.FWB.conservation.get(M.targ.chr(tidx),M.targ.start(tidx),M.targ.end(tidx)));
    conservation_track(conservation_track==200) = NaN;  % missing data
    context_and_effect_track = double(M.FWB.context_and_effect.get(M.targ.chr(tidx),M.targ.start(tidx),M.targ.end(tidx)));
    if ~strcmpi(basewise_coverage_fwb_file,'IMPUTE_FULL_COVERAGE')
      coverage_track = double(M.FWB.basewise_coverage.get(M.targ.chr(tidx),M.targ.start(tidx),M.targ.end(tidx)));
    else
      coverage_track = ones(genelength,1);
    end

    % impute at least median coverage in radius around each mutation
    medcov = median(coverage_track(coverage_track>0));
    for i=1:nm,m=midx(i);
      i1 = max(1,M.mut.trackpos(m)-P.radius_to_impute_coverage_around_mutations);
      i2 = min(genelength,M.mut.trackpos(m)+P.radius_to_impute_coverage_around_mutations);
      coverage_track(i1:i2) = max(coverage_track(i1:i2),medcov);
    end
  
    % simplify coverage track (reduce number of discrete coverage levels)
    maxcov = max(coverage_track);
    if maxcov > P.max_coverage_bins
      coverage_track_factor = maxcov / P.max_coverage_bins;
      coverage_track = round(coverage_track / coverage_track_factor); % maybe should be ceil?
      maxcov = max(coverage_track);
    end
    if maxcov==0, break; end   % this gene has no coverage

    % enumerate throwable positions for each mutation flavor
    categ = M.mut.categ_idx(midx);
    effect = M.effect.permutations_effect_idx(M.mut.effect_idx(midx));
    categ(effect==5) = inf;  % collapse all indels to categ=inf
    [uce tmp throw_flavor] = unique([categ effect],'rows');  
    nflavors = size(uce,1); flavor_counts = histc(throw_flavor,1:nflavors);
    throwable = cell(nflavors,1);
    for fi=1:nflavors
      c = uce(fi,1); e = uce(fi,2);
      if e<5          % e=3=missense or e=4=null=nonsense/splice: need to match effect in Q
        [context_and_effect_idx newbase_idx] = find(M.context_and_effect.Q(:,c,:)==e);
        cn = [context_and_effect_idx newbase_idx];
        % make sure the classes of the observed positions are included (in case track has disagreements)
        mm = midx(throw_flavor==fi);
        cn_obs = [context_and_effect_track(M.mut.trackpos(mm)) M.mut.newbase_idx(mm)];
        cn = unique([cn;cn_obs],'rows');
        % make list of throwable for this flavor
        nx = size(cn,1); x = cell(nx,1);
        for i=1:nx
          idx = find(context_and_effect_track==cn(i,1));
          y = cell(length(idx),1);
          for j=1:length(idx), y{j} = repmat([idx(j) cn(i,2)],coverage_track(idx(j)),1); end
          x{i} = cat(1,y{:});
        end
        throwable{fi} = cat(1,x{:});
      else            % e=5=indel: can throw to whole territory, weighted by the coverage track
        x = cell(genelength,1);
        for i=1:genelength, x{i}=i*ones(coverage_track(i),1); end
        x = cat(1,x{:});
        throwable{fi} = [x 5*ones(size(x))]; % newbase=5
      end
    end
    nthrowable = cellfun('length',throwable);
    if sum(nthrowable)==0, break; end
  
    % range of metrics and bins for joint distrib
    min_clust = 0; max_clust = 1; min_cons = min(conservation_track); max_cons = max(conservation_track);
    if isnan(min_cons), min_cons = 0; end
    if isnan(max_cons), max_cons = 100; end
    nbins = 100; binsize_clust = (max_clust-min_clust)/(nbins-1); binsize_cons = (max_cons-min_cons)/(nbins-1);
    
    % calculate metrics for the observed mutations
    obs_cons = nanmean(conservation_track(M.mut.trackpos(midx)));
    if isnan(obs_cons), obs_cons = 50; end
    obs_cons = obs_cons / P.permutations_min_effect_size;
    obs_cons_bin = 1+floor((obs_cons - min_cons)/binsize_cons);
    obs_clust = new_clustering_statistic(M.mut.trackpos(midx),genelength,P.clustering_metric,true) / P.permutations_min_effect_size;
    obs_clust_bin = 1+floor((obs_clust - min_clust)/binsize_clust);
  
    %%%%%%%%%%%%%%%%%%%
    % PERMUTATIONS
    %%%%%%%%%%%%%%%%%%%
    
    k_clust = 0; k_cons = 0; joint_hist = zeros(nbins,nbins);
    nperm = 0; first_check = 1000; check_every = 10000; finished = false; tt1=tic;
    while(~finished), nperm = nperm + 1;
      % randomly throw mutations
      thrown_muts = cell(nflavors,1);
      for j=1:nflavors, if nthrowable(j)>0, thrown_muts{j} = throwable{j}(ceil(nthrowable(j)*rand(flavor_counts(j),1)),:);end,end
      perm_mutpos = cat(1,thrown_muts{:});
      
      % calculate metrics for this permutation and increment histograms
      perm_cons = nanmean(conservation_track(perm_mutpos(:,1)));
      if isnan(perm_cons), perm_cons = 50; end
      perm_clust = new_clustering_statistic(perm_mutpos(:,1),genelength,P.clustering_metric);
      if perm_cons>=obs_cons, k_cons=k_cons+1; end
      if perm_clust>=obs_clust, k_clust=k_clust+1; end
      bin_cons = 1+floor((perm_cons - min_cons)/binsize_cons);
      bin_clust = 1+floor((perm_clust - min_clust)/binsize_clust);
      joint_hist(bin_clust,bin_cons) = joint_hist(bin_clust,bin_cons) + 1;
      
      if nperm~=first_check && mod(nperm,check_every)>0 && nperm<P.maxperm, continue; end  % don't need to check p-values every iteration
        
      % calculate marginal and joint p-values
      [pclust ci_ratio_clust] = calc_pval_and_ci_ratio(k_clust,nperm);
      [pcons ci_ratio_cons] = calc_pval_and_ci_ratio(k_cons,nperm);
      landfilled_hist = landfill(joint_hist/nperm);
      bin_score = -log10(landfilled_hist);
      obs_score = bin_score(obs_clust_bin,obs_cons_bin);
      k_joint = sum(joint_hist(bin_score>=obs_score));
      [p_joint ci_ratio_joint] = calc_pval_and_ci_ratio(k_joint,nperm);
      max_ci_ratio = max([ci_ratio_joint,ci_ratio_cons,ci_ratio_clust]);
      finished = (max_ci_ratio<=P.theta) | (nperm>=P.maxperm);
      
      % STATUS REPORT
      sec_remaining = toc(tt) * (M.ng-g)/g;
      eta = now + (sec_remaining/(60*60*24));
      fprintf('%s  %s %5d/%5d %-11s %-4d %7d  CV %-0.6f CL %-0.6f FN %-0.6f  [ %-4d %-4d %-4d %-4d %-4d %-4d ]\n',...
              datestr(now,'HH:MM:SS'),datestr(eta,'HH:MM:SS'),g,M.ng,M.gene.name{g},nm,nperm,M.gene.pCVmax(g),pclust,pcons,...
              qqhist(1),qqhist(2),qqhist(3),qqhist(4),qqhist(5),qqhist(6));
      gene_status_reported = true;
      
    end  % next permutation
    
    % record results of permutations
    M.gene.nperm(g) = nperm;
    M.gene.pCL(g) = pclust;
    M.gene.pFN(g) = pcons;
    M.gene.pCLFN(g) = p_joint;
    
    break
  end   % end of container loop

  if ~gene_status_reported
      % STATUS REPORT
      sec_remaining = toc(tt) * (M.ng-g)/g;
      eta = now + (sec_remaining/(60*60*24));
      fprintf('%s  %s %5d/%5d %-11s %-4d %7d  CV %-0.6f CL %-0.6f FN %-0.6f  [ %-4d %-4d %-4d %-4d %-4d %-4d ]\n',...
              datestr(now,'HH:MM:SS'),datestr(eta,'HH:MM:SS'),g,M.ng,M.gene.name{g},nm,nperm,pCVmax,pclust,pcons,...
              qqhist(1),qqhist(2),qqhist(3),qqhist(4),qqhist(5),qqhist(6));
      gene_status_reported = true;
  end

  % increment text representation of q-q plot for status report
  sp = min(6,floor(-log10(min([pCVmax pclust pcons]))));
  if sp>=1, qqhist(sp)=qqhist(sp)+1; end

 catch me
   fprintf('ERROR with gene %s\n',M.gene.name{g});
   disp(me);
   disp(me.message);
   %keyboard
 end

end, fprintf('\n');   % next gene

% close track files
M.FWB.conservation.close();
M.FWB.context_and_effect.close();
if ~strcmpi(basewise_coverage_fwb_file,'IMPUTE_FULL_COVERAGE')
  M.FWB.basewise_coverage.close();
end

% cap p-values
pcap = 10.^-(M.gene.nmut); M.gene.pFN = max(M.gene.pFN,pcap); M.gene.pCL = max(M.gene.pCL,pcap); M.gene.pCLFN = max(M.gene.pCLFN,pcap);
pcap = 1./M.gene.nperm;    M.gene.pFN = max(M.gene.pFN,pcap); M.gene.pCL = max(M.gene.pCL,pcap); M.gene.pCLFN = max(M.gene.pCLFN,pcap);
pcap = 1e-16;              M.gene.pCVmin = max(M.gene.pCVmin,pcap); M.gene.pCVmid = max(M.gene.pCVmid,pcap); M.gene.pCVmax = max(M.gene.pCVmax,pcap);

% combine p-values
M.gene.pmid = max(1e-16,fisher_combined_p([M.gene.pCVmid M.gene.pCLFN]));
M.gene.pmax = max(1e-16,fisher_combined_p([M.gene.pCVmax M.gene.pCLFN]));

% for genes that didn't undergo permutations, just take the MutSigCV values
idx = find(isnan(M.gene.pCLFN)|isnan(M.gene.nperm)|M.gene.nperm==0);
M.nperm(idx)=0;
M.gene.pFN(idx)=nan; M.gene.pCL(idx)=nan; M.gene.pCLFN(idx)=nan;
M.gene.pmid(idx) = M.gene.pCVmid(idx); M.gene.pmax(idx) = M.gene.pCVmax(idx);

% FDR
M.gene.q = calc_fdr_value(M.gene.pmax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% WRITE RESULTS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%d genes with q<=0.1\n',sum(M.gene.q<=0.1));

fprintf('Saving results... ');

% save categories file
save_struct(M.categ,[output_dir '/mutcategs.txt']);

% save patient counts and rates file
save_struct(rmfields(M.pat, {'N_c' 'n_c' 'rate_c'}),[output_dir '/patient_counts_and_rates.txt']);

% save final mutation list
% maintain MAF spec, duplicate fields.
if strcmp(build, 'hg18'),
  M.mut.pos = M.mut.pos18;
  M.mut = rmfield(M.mut, 'pos18');
end
M.mut.categ = M.mut.categ_idx;
M.mut = reorder_struct(M.mut, ~isnan(M.mut.context65));
tmp = keep_fields(M.mut, {'gene', 'chr', 'pos', 'type', 'classification', 'ref_allele', 'patient'});
M.mut = rename_fields(M.mut, ...
	{'gene', 'chr', 'pos', 'type', 'classification', 'ref_allele', 'patient'}, ...
	{'Hugo_Symbol', 'Chromosome', 'Start_Position', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Sample_Barcode'});
M.mut = merge_structs({M.mut, tmp});
save_struct(M.mut,[output_dir '/final_analysis_set.maf']);

% SAVE full output matfile
[G Go] = sort_struct(M.gene,{'pmax','npat'},[1 -1]);
G = rename_field(G,{'name','pCVmax','pmax','nsyn','nnon'},{'gene','pCV','p','nsil','nstp'});
G.nnon = G.nmis+G.nstp+G.nspl+G.nind;
G = move_field_to_before(G,'nnon','npat');
save([output_dir '/results.mat'],'G');

% SAVE sig_genes table
G.rank = as_column(1:slength(G));
g = keep_fields_if_exist(G,{'rank', 'gene','longname','codelen','nnei','nncd','nsil','nmis','nstp',...
		    'nspl','nind','nnon','npat','nsite','pCV','pCL','pFN','p','q'});
save_struct(g,[output_dir '/sig_genes.txt']);

%SAVE sample_sig_genes table
M.gene.ssg = zeros(slength(M.gene), M.np); %note this doesn't conform to 1 column vector per gene standard, but at this point that doesn't matter.
M.mut.is_nonsil = M.mut.is_missense | M.mut.is_nonsense | M.mut.is_splice | M.mut.is_indel;
for gi = 1:M.ng,
  h = histc(M.mut.pat_idx(M.geneidx{gi}).*M.mut.is_nonsil(M.geneidx{gi}), 1:M.np);
  if isempty(h), continue; end
  M.gene.ssg(gi, :) = h;
end
G.ssg = M.gene.ssg(Go, :);

%now print to file
ssgfile = fopen([output_dir '/sample_sig_gene_table.txt'], 'wt');

%header
fprintf(ssgfile, 'gene');
for i = 1:M.np,
  fprintf(ssgfile, '\t%s', M.pat.name{i});
end

fprintf(ssgfile, '\n');

%rows
for i = 1:M.ng,
  fprintf(ssgfile, '%s\t', G.gene{i});
  fprintf(ssgfile, '%d\t', G.ssg(i, :));
  
  if G.q(i) >= 0.2, break; end %old sample_sig_gene_table cutoff

  fprintf(ssgfile, '\n');
end

fclose(ssgfile);

%SAVE per_gene.mutation_counts
pgmc = keep_fields(M.gene, {'name', 'longname', 'chr', 'gene_start', 'gene_end', 'tot_exon_len', 'gc'});
pgmc = rename_field(pgmc, {'gene_start', 'gene_end', 'tot_exon_len'}, {'start', 'end', 'len'});
pgmc.cov_gidx = as_column(1:M.ng);

%header fieldnames
pgmc_fn1 = {'name', 'cov_gidx', 'longname', 'chr', 'start', 'end', 'len', 'gc'};
pgmc_fn = cat(2, pgmc_fn1, {M.pat.name{:}});

%now print to file 
pgmcfile = fopen([output_dir '/per_gene.mutation_counts.txt'], 'wt');

%header
for i = 1:length(pgmc_fn),
  fprintf(pgmcfile, '%s\t', pgmc_fn{i});
end
fprintf(pgmcfile, '\n');

%rows
for i = 1:M.ng,
  %is easier to just print each struct element separately
  for j = 1:length(pgmc_fn1),
    if iscell(pgmc.(pgmc_fn1{j})),
      fprintf(pgmcfile, '%s\t', pgmc.(pgmc_fn1{j}){i});
    else
      fprintf(pgmcfile, '%d\t', pgmc.(pgmc_fn1{j})(i));
    end
  end

  fprintf(pgmcfile, '%d\t', M.gene.ssg(i, :));
  fprintf(pgmcfile, '\n');
end

fclose(pgmcfile);

fprintf('Done.\n');
