function [pval, score, aux] = calculate_significance(N_work,n_work,mu_work,P)
% [pval, score, aux] = calculate_significance(N_work,n_work,mu_work,P)
%
% input matrices:  rows=genes, cols=categories, pages=patients
%
% if mu_work is supplied as {X_work,x_work}, will use beta-binomial distribution
%
% Mike Lawrence 2008-2011

if nargin<3, error('requires at least three arguments'); end

%if nargout<2, error('output format of calculate_significance has changed'); end

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'sig_calculation_method','projection');

if strcmpi(P.sig_calculation_method,'projection')
  P=impose_default_value(P,'null_boost_factor',1);
  P=impose_default_value(P,'projection_bin_cap',1e-6);
  P=impose_default_value(P,'use_sample_specific_mutation_rates',true);
  if ~P.use_sample_specific_mutation_rates
    error('Projection method uses sample-specific mutation rates by definition');
  end
  P=impose_default_value(P,'projection_use_new_style_scores',false);
  P=impose_default_value(P,'null_boost_factor',1);
elseif strcmpi(P.sig_calculation_method,'projection_2D')
  P=impose_default_value(P,'null_boost_factor',1);
  P=impose_default_value(P,'projection_bin_cap',1e-6);
  P=impose_default_value(P,'use_sample_specific_mutation_rates',true);
  if ~P.use_sample_specific_mutation_rates
    error('Projection method uses sample-specific mutation rates by definition');
  end
  P=impose_default_value(P,'projection_use_new_style_scores',true);
  P=impose_default_value(P,'null_boost_factor',1);
elseif strcmpi(P.sig_calculation_method,'concatenation')
  P=impose_default_value(P,'skip_landfilling',false);
  P=impose_default_value(P,'use_sample_specific_mutation_rates',false);
  P=impose_default_value(P,'convolution_maxmuts',1000);
  P=impose_default_value(P,'convolution_maxmuts_limit',100000);
  P=impose_default_value(P,'null_boost_factor',1);
elseif strcmpi(P.sig_calculation_method,'concatenation2')
  P=impose_default_value(P,'convolution_maxmuts',1000);
  P=impose_default_value(P,'convolution_maxmuts_limit',100000);
  P=impose_default_value(P,'null_boost_factor',1);
elseif strcmpi(P.sig_calculation_method,'concatenation_feb2012')
  P=impose_default_value(P,'use_sample_specific_mutation_rates',false);
  P=impose_default_value(P,'skip_landfilling',false);
  P=impose_default_value(P,'convolution_maxmuts',1000);
  P=impose_default_value(P,'convolution_maxmuts_limit',100000);
  P=impose_default_value(P,'semi_exact_procedure_bin_size',1);
elseif strcmpi(P.sig_calculation_method,'concatenation_1D')
  P=impose_default_value(P,'use_sample_specific_mutation_rates',false);
  %(ok)
elseif strcmpi(P.sig_calculation_method,'bestcat')
  P=impose_default_value(P,'use_sample_specific_mutation_rates',true);
  P=impose_default_value(P,'convolution_maxmuts',1000);
  P=impose_default_value(P,'convolution_maxmuts_limit',100000);
else
  error('Unsupported significance calculation method "%s"', P.sig_calculation_method);  
end

P=impose_default_value(P,'convolution_minbins',10000);
P=impose_default_value(P,'min_score_to_proceed_to_convolutions',-inf);
P=impose_default_value(P,'mutsig_projection_tolerance_cutoff',0);

if isfield(P,'semi_exact_procedure_bin_size')
  fprintf('NOTE: parameter no longer used: semi_exact_procedure_bin_size\n');
end

P=impose_default_value(P,'simulation',false);

oldflds = {'use_1D_concatenation','use_semi_exact_procedure','mutation_rate_to_use','manual_mutation_rate',...
  'zero_top_genes','top_genes_to_zero','patient_subset'};
for i=1:length(oldflds), if isfield(P,oldflds{i}), error('No longer supported: P.%s',oldflds{i}); end, end

if ndims(N_work)~=ndims(n_work), error('input matrices should be of like dimensionality'); end
if any(size(N_work)~=size(n_work)), error('inconsistent input matrix sizes'); end
[ng ncat np] = size(N_work);

if isfield(P,'first_gene_to_calculate') && isfield(P,'last_gene_to_calculate')
  gta = P.first_gene_to_calculate:P.last_gene_to_calculate;
elseif isfield(P,'first_gene_to_calculate') || isfield(P,'last_gene_to_calculate')
  error('Must set both or neither of P.first_gene_to_calculate and P.last_gene_to_calculate');
elseif isfield(P,'genes_to_analyze') && ~isempty(P.genes_to_analyze)
  gta = P.genes_to_analyze;
  if ischar(gta), gta={gta}; end
  if islogical(gta), gta=find(gta); end
  if iscellstr(gta)
    if isfield(P,'gene_names')
      gta = listmap(gta,P.gene_names);
      gta(isnan(gta)) = [];
    else
      error('need P.gene_names if P.genes_to_analyze lists gene names');
    end
  elseif isnumeric(gta)
    % OK
  else
    error('unknown format for P.genes_to_analyze');
  end
else
  gta = 1:ng;
end
gta = as_row(gta);

if isfield(P,'gene_name') && ~isfield(P,'gene_names'), P = rename_field(P,'gene_name','gene_names'); end
if isfield(P,'patient_name') && ~isfield(P,'patient_names'), P = rename_field(P,'patient_name','patient_names'); end
P=impose_default_value(P,'gene_names',str2cell(sprintf('gene%d\n',1:ng)));
P=impose_default_value(P,'patient_names',str2cell(sprintf('patient%d\n',1:np)));

% is mu_work supplied as {X_work,x_work} for beta-binomial?
if isnumeric(mu_work)
  using_betabinom = false;
  if ndims(mu_work)~=ndims(N_work), error('input matrices should be of like dimensionality'); end
  if any(size(N_work)~=size(mu_work)), error('inconsistent input matrix sizes'); end
  % check to make sure we won't need to take log(0)
  [a b c] = find(n_work>0 & (N_work==0 | mu_work==0));
  if ~isempty(a)
    fprintf('Zero mutrate and/or coverage in %d bin(s) with nonzero mutation counts: removing these impossible mutations.\n',length(a));
    n_work(N_work==0 | mu_work==0) = 0;
  end
elseif iscell(mu_work) && length(mu_work)==2
  using_betabinom = true;
  X_work = mu_work{1};
  x_work = mu_work{2};
  clear mu_work;
  if ndims(X_work)~=ndims(N_work) ||ndims(X_work)~=ndims(x_work), error('input matrices should be of like dimensionality'); end
  if any(size(N_work)~=size(x_work)) || any(size(N_work)~=size(X_work)), error('inconsistent input matrix sizes'); end
  if any(x_work(:)>X_work(:)), error('x_work and X_work reversed?'); end
  fprintf('using beta-binomial distribution\n');
  % check to make sure we won't need to take log(0)
  [a b c] = find(n_work>0 & N_work==0);
  if ~isempty(a)
    fprintf('Zero coverage in %d bin(s) with nonzero mutation counts: removing these impossible mutations.\n',length(a));
    n_work(N_work==0) = 0;
  end
else
  error('unknown format for mu_work');
end

%% make sure n,N,X are integers so we don't have problems with binopdf
n_work = round(n_work);
N_work = round(N_work);
if using_betabinom
  X_work = round(X_work);
end

%%%% SIMULATOR
%%%% replace observed data with simulated data?
if P.simulation
  if using_betabinom
    n_work = replace_with_simulated_data(N_work,{X_work,x_work},P);
  else
    n_work = replace_with_simulated_data(N_work,mu_work,P);
  end
end

%%% SIGNIFICANCE CALCULATION
%%% which method to use?
if strcmpi(P.sig_calculation_method,'concatenation_1D')
  METHOD = 0;
elseif strcmpi(P.sig_calculation_method, 'concatenation')
  if ~P.use_sample_specific_mutation_rates
    METHOD = 1;
  else % sample-specific mutation rates
    METHOD = 2;
  end
elseif strcmpi(P.sig_calculation_method, 'concatenation2')
  if P.use_sample_specific_mutation_rates
    METHOD = 0.98;
  else % no sample-specific mutation rates
    error('NOT SUPPORTED');
  end
elseif strcmpi(P.sig_calculation_method, 'concatenation_feb2012')
  if ~P.use_sample_specific_mutation_rates
    METHOD = 0.95;
  else % sample-specific mutation rates
    error('NOT SUPPORTED');
  end
elseif strcmpi(P.sig_calculation_method, 'projection')
  METHOD = 3;
elseif strcmpi(P.sig_calculation_method, 'projection_2D')
  METHOD = 3.5;
elseif strcmpi(P.sig_calculation_method, 'bestcat')
  METHOD = 4;
else
  error('Unsupported significance calculation method "%s"', P.sig_calculation_method);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Calculating significance of each gene.\n');

pval = nan(ng,1);
score = nan(ng,1);
aux = [];

switch METHOD

case 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1D concatenation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ncat>1 || np>1
  error('for 1D concatenation, please collapse categories and patients before calling calculate_significance');
end

fprintf('Using 1D concatenation.\n');

if using_betabinom

  pval = nan(length(n_work),1);
  for i=1:length(n_work), if ~mod(i,1000), fprintf('%d/%d ',i,length(n_work)); end
    pval(i) = 1-hyge2cdf(n_work(i)-1,N_work(i),x_work(i),X_work(i));
  end, fprintf('\n');

else

  pval = 1-binocdf(n_work-1,N_work,mu_work);

end

score = -log10(pval);

case 0.95
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% landfilled concatenation-with-categories method
% feb2012 freeze
%
% using global BMR
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('Using landfilled concatenation-with-categories method with global BMR.\n');

% collapse across samples

N_work = sum(N_work,3);
n_work = sum(n_work,3);

% verify that mu_work is the same across patients
if np>1
  all_same = true;
  for i=2:np
    if using_betabinom
      flag = any(x_work(:,:,i)~=x_work(:,:,1) | ...
                 X_work(:,:,i)~=X_work(:,:,1));
    else
      flag = any(mu_work(:,:,i)~=mu_work(:,:,1));
    end
    if flag
      fprintf('WARNING:  Not all patients have same BMRs\n');
      all_same = false;
      break;
    end
  end
  if all_same
    if using_betabinom
      x_work = x_work(:,:,1);
      X_work = X_work(:,:,1);
    else
      mu_work = mu_work(:,:,1);
    end
  else
    fprintf('Taking mean across samples.\n');
    if using_betabinom
      x_work = mean(x_work,3);
      X_work = mean(X_work,3);
    else
      mu_work = mean(mu_work,3);
    end
  end
end

fprintf('gene ');

for g=gta
  if mod(g,1000)==0, fprintf('%d/%d ', g, ng); end
  gname = P.gene_names{g};

  % save time if gene has no mutations
  if fullsum(n_work(g,1:ncat))==0
    pval(g) = 1;
    continue
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 1
  % prepare table of piecewise probabilities and scores
  % by mutation category and number of mutations
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  MAXMUTS = 200;
  if any(n_work(g,:)*2>MAXMUTS)
    old = MAXMUTS;
    MAXMUTS = max(n_work(g,:))*2;
    fprintf('Increasing MAXMUTS from %d to %d to cover number of observed mutations in %s\n',old,MAXMUTS,gname);
  end

  assign_p_zero = false;
  while(true)   % (may need to re-try with increased MAXMUTS)
    if MAXMUTS > P.convolution_maxmuts_limit
      fprintf('P.convolution_maxmuts_limit exceeded: assigning p=0\n');
      assign_p_zero = true;
      break;
    end

    prob_dist = zeros(ncat,MAXMUTS+1);    % second index is n+1 (to accommodate n=0)
    score_dist = zeros(ncat,MAXMUTS+1);
 
    need_to_increase = false;
    for c=1:ncat
      N = N_work(g,c);

      if using_betabinom
        x = x_work(g,c);
        X = X_work(g,c);
        nexp = floor(N*x/X);
        prob_dist(c,:) = hyge2pdf(0:MAXMUTS,N,x,X);
      else
        mu = mu_work(g,c);
        nexp = floor(N*mu);
        %      prob_dist(c,:) = binopdf(0:MAXMUTS,N,mu);
        prob_dist(c,:) = poisspdf(0:MAXMUTS,N*mu);
      end

      score_dist(c,:) = -log10(prob_dist(c,:));
      if ~P.skip_landfilling
        score_dist(c,1) = 0;   % no score for zero mutations
        if nexp>MAXMUTS
          need_to_increase = true;
          break;
        end
        if nexp>=1   % enforce one-sided scores ("landfill")        
          score_dist(c,1:nexp) = score_dist(c,nexp+1);
        end
      end
    end

    if need_to_increase
      old = MAXMUTS;
      MAXMUTS = nexp * 2;
      fprintf('Increasing MAXMUTS from %d to %d to cover number of observed mutations in %s\n',old,MAXMUTS,gname);
      continue;
    end

%      max_score = 500;  % <----- may need to change this treatment: can cause problem with summed scores
%      score_dist(score_dist>max_score) = max_score;
  % NOW WE JUST LEAVE IT AS Inf, and everything works out OK

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 2
    % calculate score of the observation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    score_obs = 0;
    for c=1:ncat
      n  = n_work(g,c);
      score_obs = score_obs + score_dist(c,n+1);
    end

    if isinf(score_obs)
      fprintf('infinite score: assigning p=0\n');
      assign_p_zero = true;
      break;
    end

    % compute nmax for each category:
    %  the number of mutations where score_dist meets/exceeds score_obs
    need_to_retry = false;
    nmax = nan(ncat,1);
    for c=1:ncat
      idx = find(score_dist(c,:)>=score_obs,1);
      if isempty(idx)
        old = MAXMUTS;
        MAXMUTS = MAXMUTS * 2;
        fprintf('Increasing MAXMUTS from %d to %d to cover score_obs of %s\n',old,MAXMUTS,gname);
        need_to_retry = true;
        break;
      end
      nmax(c) = idx-1;
    end
    maxnmax = max(nmax);
    
    if ~need_to_retry, break; end
  end

  if assign_p_zero
    pval(g) = 0;
    score(g) = inf;
    continue; % next gene
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 3
  % sequential convolution of mutation categories
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% OLD STYLE 

  % set up histogram
  % first: adjust binsize to make sure there is an integral number of bins per score_obs
  binsize = P.semi_exact_procedure_bin_size;
  binsize = score_obs / ceil(score_obs / binsize);
  numbins = ceil(score_obs / binsize);   % bin1 = zero only

%% NEW STYLE

%  minbins = P.convolution_minbins;
%  if score_obs > minbins
%    binsize = 1;
%  else % even for very low-scoring genes, still make sure we use at least minbins bins
%    binsize = score_obs / minbins;
%  end
%  numbins = ceil(score_obs / binsize);
%  binsize = score_obs / numbins;





  if isnan(numbins) || isinf(numbins)
    fprintf('ERROR with numbins, while processing gene %s\n', gname);
    keyboard
  end
  H = zeros(numbins,1);
  H(1) = 1;    % initial condition: all probability is in first bin (score=0, P=1)

  % sequential convolution
  offset = min(numbins,round(score_dist/binsize));
  newH = zeros(numbins,maxnmax+1);
  for c=1:ncat
    newH(:) = 0;
    for n=0:nmax(c)
      o = offset(c,n+1);
      newH(o+1:end,n+1) = prob_dist(c,n+1) .* H(1:end-o);
    end
    H = sum(newH(:,1:nmax(c)+1),2);
  end

  Pbulk = sum(H);
  pval(g) = 1 - Pbulk;
  score(g) = score_obs;

 end   % next gene

case 0.98
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% concatenation-with-categories VERSION 2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Using concatenation-with-categories VERSION 2.\n');

silence = 0;
debug = 0;

fprintf('gene ');

for g=gta

  if mod(g,1000)==0, fprintf('%d/%d ', g, ng); end

  gname = P.gene_names{g};

  % save time if gene has no mutations
  if fullsum(n_work(g,1:ncat,:))==0
    pval(g) = 1;
    continue
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 1
  % prepare table of piecewise probabilities and scores
  % by mutation category and number of mutations
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  MAXMUTS = 100;
  if any(any(n_work(g,:,:)*2>MAXMUTS))
    old = MAXMUTS;
    MAXMUTS = max(max(n_work(g,:,:)))*2;
    fprintf('Increasing MAXMUTS from %d to %d to cover number of observed mutations in %s\n',old,MAXMUTS,gname);
  end

  assign_p_zero = false;
  while(true)   % (may need to re-try with increased MAXMUTS)
    if MAXMUTS > P.convolution_maxmuts_limit
      fprintf('P.convolution_maxmuts_limit exceeded: assigning p=0\n');
      assign_p_zero = true;
      break;
    end

    prob_dist = zeros(np,ncat,MAXMUTS+1);    % second index is n+1 (to accommodate n=0)
    score_dist = zeros(np,ncat,MAXMUTS+1);

    need_to_increase = false;
    for p=1:np, for c=1:ncat
      N = N_work(g,c,p);

      if using_betabinom
        error('not supported yet');
      else
        mu = mu_work(g,c,p);
        nexp = floor(N*mu);
        score_dist(p,c,:) = (0:MAXMUTS).*-log10(mu) + (N-(0:MAXMUTS)).*-log10(1-mu); 
        prob_dist(p,c,:) = binopdf(0:MAXMUTS,N,mu);
      end

      if nexp>MAXMUTS
        need_to_increase = true;
        break;
      end
    end, end  % next c, next p

    if need_to_increase
      old = MAXMUTS;
      MAXMUTS = nexp * 2;
      fprintf('Increasing MAXMUTS from %d to %d to cover number of observed mutations in %s\n',old,MAXMUTS,gname);
      continue;
    end

%    keyboard
    score_dist = bsxfun(@minus,score_dist,score_dist(:,:,1));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 2
    % calculate score of the observation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    score_obs = 0;
    for p=1:np, for c=1:ncat
      n  = n_work(g,c,p);
      score_obs = score_obs + score_dist(p,c,n+1);
    end, end

    if isinf(score_obs)
     fprintf('infinite score: assigning p=0\n');
      assign_p_zero = true;
      break;
    end

    % compute nmax for each category:
    %  the number of mutations where score_dist meets/exceeds score_obs
    need_to_retry = false;
    nmax = nan(np,ncat);
    for p=1:np, for c=1:ncat
      idx = find(score_dist(p,c,:)>=score_obs,1);
      if isempty(idx)
        old = MAXMUTS;
        MAXMUTS = MAXMUTS * 2;
        fprintf('Increasing MAXMUTS from %d to %d to cover score_obs of %s\n',old,MAXMUTS,gname);
        need_to_retry = true;
        break;
      end
      nmax(p,c) = idx-1;
    end, end
    maxnmax = max(max(nmax));

    if ~need_to_retry, break; end
  end

  if assign_p_zero
    pval(g) = 0;
    score(g) = inf;
    continue; % next gene
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 3
  % sequential convolution of mutation categories
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % set up histogram
  % first: adjust binsize to make sure there is an integral number of bins per score_obs
  minbins = P.convolution_minbins;
  disp(minbins)
  if score_obs > minbins
    binsize = 1;
  else % even for very low-scoring genes, still make sure we use at least minbins bins
    binsize = score_obs / minbins;
  end
  numbins = ceil(score_obs / binsize);
  binsize = score_obs / (numbins-0.5);
  if isinf(numbins)||isnan(numbins), fprintf('what?\n'); keyboard;end

  H = zeros(numbins,1);
  H(1) = 1;    % initial condition: all probability is in first bin (score=0, P=1)

  % sequential convolution
  offset = min(numbins,round(score_dist/binsize));
  newH = zeros(numbins,maxnmax+1);
  for p=1:np, for c=1:ncat
    newH(:) = 0;
    for n=0:nmax(p,c)
      o = offset(p,c,n+1);
      newH(o+1:end,n+1) = prob_dist(p,c,n+1) .* H(1:end-o);
    end
    H = sum(newH(:,1:nmax(p,c)+1),2);
 end, end

  Pbulk = sum(H);
  pval(g) = 1 - Pbulk;
  score(g) = score_obs;

end   % next gene
%%%% end of method


case 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% landfilled concatenation-with-categories method 2010/12/19
%
% using global BMR
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Using landfilled concatenation-with-categories method with global BMR.\n');

% collapse across samples

N_work = sum(N_work,3);
n_work = sum(n_work,3);

% verify that mu_work is the same across patients
if np>1
  all_same = true;
  for i=2:np
    if using_betabinom
      flag = any(x_work(:,:,i)~=x_work(:,:,1) | ...
                 X_work(:,:,i)~=X_work(:,:,1));
    else
      flag = any(mu_work(:,:,i)~=mu_work(:,:,1));
    end
    if flag
      fprintf('WARNING:  Not all patients have same BMRs\n');
      all_same = false;
      break;
    end
  end
  if all_same
    if using_betabinom
      x_work = x_work(:,:,1);
      X_work = X_work(:,:,1);
    else
      mu_work = mu_work(:,:,1);
    end
  else
    fprintf('Taking mean across samples.\n');
    if using_betabinom
      x_work = mean(x_work,3);
      X_work = mean(X_work,3);
    else
      mu_work = mean(mu_work,3);
    end
  end
end

fprintf('gene ');

for g=gta
  if mod(g,1000)==0, fprintf('%d/%d ', g, ng); end
  gname = P.gene_names{g};

  % save time if gene has no mutations
  if fullsum(n_work(g,1:ncat))==0
    pval(g) = 1;
    continue
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 1
  % prepare table of piecewise probabilities and scores
  % by mutation category and number of mutations
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  MAXMUTS = 200;
  if any(n_work(g,:)*2>MAXMUTS)
    old = MAXMUTS;
    MAXMUTS = max(n_work(g,:))*2;
    fprintf('Increasing MAXMUTS from %d to %d to cover number of observed mutations in %s\n',old,MAXMUTS,gname);
  end

  assign_p_zero = false;
  while(true)   % (may need to re-try with increased MAXMUTS)
    if MAXMUTS > P.convolution_maxmuts_limit
      fprintf('P.convolution_maxmuts_limit exceeded: assigning p=0\n');
      assign_p_zero = true;
      break;
    end

    prob_dist = zeros(ncat,MAXMUTS+1);    % second index is n+1 (to accommodate n=0)
    score_dist = zeros(ncat,MAXMUTS+1);
 
    need_to_increase = false;
    for c=1:ncat
      N = N_work(g,c);

      if using_betabinom
        x = x_work(g,c);
        X = X_work(g,c);
        nexp = floor(N*x/X);
        prob_dist(c,:) = hyge2pdf(0:MAXMUTS,N,x,X);
      else
        mu = mu_work(g,c);
        nexp = floor(N*mu);
        %      prob_dist(c,:) = binopdf(0:MAXMUTS,N,mu);
        prob_dist(c,:) = poisspdf(0:MAXMUTS,N*mu);
      end

      score_dist(c,:) = -log10(prob_dist(c,:));

      if nexp>MAXMUTS
        need_to_increase = true;
        break;
      end
      if ~P.skip_landfilling
        score_dist(c,1) = 0;   % no score for zero mutations
        if nexp>=1   % enforce one-sided scores ("landfill")        
          score_dist(c,1:nexp) = score_dist(c,nexp+1);
        end
      end
      
    end

    if need_to_increase
      old = MAXMUTS;
      MAXMUTS = nexp * 2;
      fprintf('Increasing MAXMUTS from %d to %d to cover number of observed mutations in %s\n',old,MAXMUTS,gname);
      continue;
    end

%      max_score = 500;  % <----- may need to change this treatment: can cause problem with summed scores
%      score_dist(score_dist>max_score) = max_score;
  % NOW WE JUST LEAVE IT AS Inf, and everything works out OK

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 2
    % calculate score of the observation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    score_obs = 0;
    for c=1:ncat
      n  = n_work(g,c);
      score_obs = score_obs + score_dist(c,n+1);
    end

    if isinf(score_obs)
      fprintf('infinite score: assigning p=0\n');
      assign_p_zero = true;
      break;
    end

    % compute nmax for each category:
    %  the number of mutations where score_dist meets/exceeds score_obs
    need_to_retry = false;
    nmax = nan(ncat,1);
    for c=1:ncat
      idx = find(score_dist(c,:)>=score_obs,1);
      if isempty(idx)
        old = MAXMUTS;
        MAXMUTS = MAXMUTS * 2;
        fprintf('Increasing MAXMUTS from %d to %d to cover score_obs of %s\n',old,MAXMUTS,gname);
        need_to_retry = true;
        break;
      end
      nmax(c) = idx-1;
    end
    maxnmax = max(nmax);
    
    if ~need_to_retry, break; end
  end

  if assign_p_zero
    pval(g) = 0;
    score(g) = inf;
    continue; % next gene
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 3
  % sequential convolution of mutation categories
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % set up histogram
  % first: adjust binsize to make sure there is an integral number of bins per score_obs
  minbins = P.convolution_minbins;
  if score_obs > minbins
    binsize = 1;
  else % even for very low-scoring genes, still make sure we use at least minbins bins
    binsize = score_obs / minbins;
  end
  numbins = ceil(score_obs / binsize);
  binsize = score_obs / numbins;
  if isinf(numbins)||isnan(numbins), fprintf('what?\n'); keyboard; continue; end

  H = zeros(numbins,1);
  H(1) = 1;    % initial condition: all probability is in first bin (score=0, P=1)

  % sequential convolution
  offset = min(numbins,round(score_dist/binsize));
  newH = zeros(numbins,maxnmax+1);
  for c=1:ncat
    newH(:) = 0;
    for n=0:nmax(c)
      o = offset(c,n+1);
      newH(o+1:end,n+1) = prob_dist(c,n+1) .* H(1:end-o);
    end
    H = sum(newH(:,1:nmax(c)+1),2);
  end

  Pbulk = sum(H);
  pval(g) = 1 - Pbulk;
  score(g) = score_obs;

 end   % next gene


case 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% landfilled concatenation-with-categories method 2010/12/19
%
% using sample-specific BMRs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Using landfilled concatenation-with-categories method with sample-specific BMRs.\n');

silence = 0;
debug = 0;

fprintf('gene ');

for g=gta

  if mod(g,1000)==0, fprintf('%d/%d ', g, ng); end

  gname = P.gene_names{g};

  % save time if gene has no mutations
  if fullsum(n_work(g,1:ncat,:))==0
    pval(g) = 1;
    continue
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 1
  % prepare table of piecewise probabilities and scores
  % by mutation category and number of mutations
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  MAXMUTS = 100;
  if any(any(n_work(g,:,:)*2>MAXMUTS))
    old = MAXMUTS;
    MAXMUTS = max(max(n_work(g,:,:)))*2;
    fprintf('Increasing MAXMUTS from %d to %d to cover number of observed mutations in %s\n',old,MAXMUTS,gname);
  end

  assign_p_zero = false;
  while(true)   % (may need to re-try with increased MAXMUTS)
    if MAXMUTS > P.convolution_maxmuts_limit
      fprintf('P.convolution_maxmuts_limit exceeded: assigning p=0\n');
      assign_p_zero = true;
      break;
    end

    prob_dist = zeros(np,ncat,MAXMUTS+1);    % second index is n+1 (to accommodate n=0)
    score_dist = zeros(np,ncat,MAXMUTS+1);

    need_to_increase = false;
    for p=1:np, for c=1:ncat
      N = N_work(g,c,p);

      if using_betabinom
        x = x_work(g,c,p);
        X = X_work(g,c,p);
        nexp = floor(N*x/X);
        prob_dist(p,c,:) = hyge2pdf(0:MAXMUTS,N,x,X);
      else
        mu = mu_work(g,c,p);
        nexp = floor(N*mu);
        %      prob_dist(c,:) = binopdf(0:MAXMUTS,N,mu);
        prob_dist(p,c,:) = poisspdf(0:MAXMUTS,N*mu);
      end

      score_dist(p,c,:) = -log10(prob_dist(p,c,:));
      score_dist(p,c,1) = 0;   % no score for zero mutations
      if nexp>MAXMUTS
        need_to_increase = true;
        break;
      end
      if nexp>=1   % enforce one-sided scores ("landfill")
        score_dist(p,c,1:nexp) = score_dist(p,c,nexp+1);
      end
    end, end  % next c, next p

    if need_to_increase
      old = MAXMUTS;
      MAXMUTS = nexp * 2;
      fprintf('Increasing MAXMUTS from %d to %d to cover number of observed mutations in %s\n',old,MAXMUTS,gname);
      continue;
    end

%      max_score = 500;  % <----- may need to change this treatment: can cause problem with summed scores
%      score_dist(score_dist>max_score) = max_score;
  % NOW WE JUST LEAVE IT AS Inf, and everything works out OK

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 2
    % calculate score of the observation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    score_obs = 0;
    for p=1:np, for c=1:ncat
      n  = n_work(g,c,p);
      score_obs = score_obs + score_dist(p,c,n+1);
    end, end

    if isinf(score_obs)
      fprintf('infinite score: assigning p=0\n');
      assign_p_zero = true;
      break;
    end

    % compute nmax for each category:
    %  the number of mutations where score_dist meets/exceeds score_obs
    need_to_retry = false;
    nmax = nan(np,ncat);
    for p=1:np, for c=1:ncat
      idx = find(score_dist(p,c,:)>=score_obs,1);
      if isempty(idx)
        old = MAXMUTS;
        MAXMUTS = MAXMUTS * 2;
        fprintf('Increasing MAXMUTS from %d to %d to cover score_obs of %s\n',old,MAXMUTS,gname);
        need_to_retry = true;
        break;
      end
      nmax(p,c) = idx-1;
    end, end
    maxnmax = max(max(nmax));

    if ~need_to_retry, break; end
  end

  if assign_p_zero
    pval(g) = 0;
    score(g) = inf;
    continue; % next gene
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 3
  % sequential convolution of mutation categories
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % set up histogram
  % first: adjust binsize to make sure there is an integral number of bins per score_obs
  minbins = P.convolution_minbins;
  disp(minbins)
  if score_obs > minbins
    binsize = 1;
  else % even for very low-scoring genes, still make sure we use at least minbins bins
    binsize = score_obs / minbins;
  end
  numbins = ceil(score_obs / binsize);
  binsize = score_obs / numbins;
  if isinf(numbins)||isnan(numbins), fprintf('what?\n'); keyboard; continue;end

  H = zeros(numbins,1);
  H(1) = 1;    % initial condition: all probability is in first bin (score=0, P=1)

%keyboard

  % sequential convolution
  offset = min(numbins,round(score_dist/binsize));
  newH = zeros(numbins,maxnmax+1);
  for p=1:np, for c=1:ncat
    newH(:) = 0;
    for n=0:nmax(p,c)
      o = offset(p,c,n+1);
      newH(o+1:end,n+1) = prob_dist(p,c,n+1) .* H(1:end-o);
    end
    H = sum(newH(:,1:nmax(p,c)+1),2);
%    fprintf('p %d    c %d    1-sum(H) %d\n',p,c,1-sum(H));
%    bar(H);
%    keyboard
 end, end

  Pbulk = sum(H);
  pval(g) = 1 - Pbulk;
  score(g) = score_obs;

end   % next gene
%%%% end of sample-specific-BMR concatenation method


case 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  PROJECTION METHOD
%
%     Analyzes mutational data taking into account the possibility that
%     multiple mutations in the same gene in the same sample may not have
%     been independent events.   -- added 2008-04-23
%
%  1. Projects mutational profile of each sample to a restricted subspace,
%     where at most one mutation is represented,
%     selected from the most unlikely category.
%  2. For each sample, calculates the probability of that sample projecting
%     to each point in the collapsed projection.
%  3. Calculates P value for the gene based on the probability of obtaining
%     a distribution of samples in the collapsed space at least as extreme
%     as the observed one.
%
%  Extensively reworked 2011-Aug/Sep
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('Using projection method.\n');
if P.null_boost_factor>1
  fprintf('\tUsing null boost.\n');
end

silence = false;
debug = true;

aux = cell(1,1);
aux{1} = nan(ng,1);  % will store #seconds elapsed per gene

%for g = find(strcmp('TTN',P.gene_names))
for g=gta

  tt = tic;

  if ~mod(g,1000), fprintf('%d/%d ', g, ng); end
  gname = P.gene_names{g};

  if debug, fprintf('gene%d: ',g); end

  % STEP 1
  % for each sample, prioritize mutation categories according to how likely
  % it would be for this gene x sample to have a mutation there by chance.

  N = reshape(N_work(g,:,:),ncat,np)';
  n = reshape(n_work(g,:,:),ncat,np)';
  N(n>N)=n(n>N);  % make sure we don't have N>n

  if using_betabinom
    x = reshape(x_work(g,:,:),ncat,np)';
    X = reshape(X_work(g,:,:),ncat,np)';
    Pmut = 1-hyge2pdf_at_zero(N,x,X);
  else
    mu = reshape(mu_work(g,:,:),ncat,np)';
    Pmut = 1-((1-mu).^N);
  end

  if any(Pmut(:)<0 | Pmut(:)>1)
    fprintf('WARNING:  Pmut<0 or Pmut>1 at some entries\n');
    keyboard;
  end

  % determine each patient's priority order of categories
  [Pmut priority] = sort(Pmut,2,'descend');
  Pclear = 1-Pmut;

  % STEP 2
  % for each sample, compute probability that it would been of each degree.

  Pdeg = nan(np,ncat+1);
  for d=0:ncat
    % has to be clear in all categories > d
    Pdeg(:,d+1) = prod(Pclear(:,d+1:end),2);
    % and nonclear in category d (if d>0)
    if d>0, Pdeg(:,d+1) = Pdeg(:,d+1) .* Pmut(:,d); end
  end

  %% STEP 2a: calculate score for a sample being of each possible degree
  if P.projection_use_new_style_scores
    Sdeg = [zeros(np,1) -log10(Pmut)];  % NEW STYLE:  score = probability of having a mutation in that category
  else
    Sdeg = -log10(Pdeg);                % OLD STYLE:  score = probability of being that degree
  end

  %% STEP 2b: cap scores to fix problem of mutations in sparsely covered bins
  PROBCAP = P.projection_bin_cap;
  scorecap = -log10(PROBCAP);
  [pi ci] = find(Sdeg>scorecap);
  if ~isempty(pi)
    flag = false;
    for i=1:length(pi)
      if ci(i)>1 && n(pi(i),ci(i)-1)>0, flag=true; end
    end
    if flag
      fprintf('In gene %s the following bins exceed PROBCAP:\n',gname);
      for i=1:length(pi)
        if ci(i)>1 && n(pi(i),ci(i)-1)>=1
          fprintf('    %s\tcateg %d',P.patient_names{pi(i)}, ci(i));
          fprintf('\tN = %-6.0f  n = %.0f',N(pi(i),ci(i)-1),n(pi(i),ci(i)-1));
          if ~using_betabinom
            fprintf('\tmu = %10.2d', mu(pi(i),ci(i)-1));
          else
            fprintf('\tX = %-6.0f  x = %.2f', X(pi(i),ci(i)-1),x(pi(i),ci(i)-1));
          end
          fprintf('\n');
        end
      end
    end
    Sdeg = min(Sdeg,scorecap);
  end
  
  % STEP 2c: landfill
  for pi=1:np
    for di=ncat-1:-1:0
      if Sdeg(pi,di+1)>=Sdeg(pi,di+2)
        Sdeg(pi,di+1) = Sdeg(pi,di+2) - eps;
  end,end,end
  Sdeg(:,1) = 0;  % no score for degree0

  % STEP 2d: add "null boost" if applicable
  if P.null_boost_factor>1
    null_boost_score = log10(P.null_boost_factor);
    nullcat = max(priority(:));  % assume null is the LAST category
    priority2 = [zeros(np,1) priority];
    Sdeg(priority2==nullcat) = Sdeg(priority2==nullcat) + null_boost_score;
  end

  % STEP 3
  % determine actual degree of each sample, and score_obs for the gene

  S = Sdeg;  % (for debugging display purposes)

  degree = zeros(np,1);
  score_obs = 0;
  for p = 1:np
    for d = ncat:-1:1
      c = priority(p,d);
      if n(p,c)>0, degree(p) = d; break; end
    end
    score_obs = score_obs + Sdeg(p,degree(p)+1);
    S(p,degree(p)+1) = -S(p,degree(p)+1);   % (for debugging display purposes)
  end

  if score_obs<0
    fprintf('What?? score_obs<0???\n');
    score(g) = 0;
    pval(g) = 1;
  elseif isnan(score_obs)
    fprintf('What?? score_obs is nan?\n');
    score(g) = 0;
    pval(g) = 1;
  elseif score_obs==0
    score(g) = 0;
    pval(g) = 1;
  elseif isinf(score_obs)
    pval(g) = 0;
    score(g) = inf;
  else

    if score_obs<P.min_score_to_proceed_to_convolutions
      score(g) = score_obs;
      pval(g) = 0.5;
      continue
    end
    
    % STEP 4
    % compute P value for gene by convolutions

    % set up histogram
    %    -- choose binsize to be an even division of score_obs
    %       so that all probability in the last bin is *less* than the observed score

    minbins = P.convolution_minbins;
    if score_obs > minbins
      binsize = 1;
    else % even for very low-scoring genes, still make sure we use at least minbins bins
      binsize = score_obs / minbins;
    end
    numbins = ceil(score_obs / binsize);
    binsize = score_obs / numbins;
    if isinf(numbins)||isnan(numbins), fprintf('what?\n'); keyboard; continue;end

    H = zeros(numbins,1);
    H(1) = 1;  % initial condition: all probability is in first bin
    % sequential convolution
    offset = min(numbins, round(Sdeg/binsize));
    newH = zeros(numbins,ncat+1);
    for p = 1:np
      newH(:) = 0;
      for d = 0:ncat
        o = offset(p,d+1);
        newH(o+1:end,d+1) = Pdeg(p,d+1) .* H(1:end-o);
      end
      H = sum(newH,2);
    end

    pval(g) = 1-sum(H);
    if pval(g)<0, pval(g) = 0; end
    score(g) = score_obs;
    aux{1}(g) = toc(tt);
  end

%  if strcmp(gname,'TP53'), keyboard;end

  if ~silence
    if debug || pval(g)<1e-4 || length(gta)<100
      fprintf('%d\t%s\t%d\n', g, gname, pval(g));
    end
  end

end % next gene

%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of projection method

case 3.5
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  PROJECTION METHOD_2D
  %
  %  takes at most two mutations from each sample
  %
  %  2012-08-29
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  fprintf('Using projection_2D method.\n');
  if P.null_boost_factor>1
    fprintf('\tUsing null boost.\n');
  end

  silence = false;
  debug = false;

  for g=gta
    if ~mod(g,1000), fprintf('%d/%d ', g, ng); end
    gname = P.gene_names{g};
    
    if debug, fprintf('gene%d: ',g); end

    % STEP 1
    % for each sample, prioritize mutation categories according to how likely
    % it would be for this gene x sample to have a mutation there by chance.
    
    N = reshape(N_work(g,:,:),ncat,np)';
    n = reshape(n_work(g,:,:),ncat,np)';
    N(n>N)=n(n>N);  % make sure we don't have N>n

    if using_betabinom
      x = reshape(x_work(g,:,:),ncat,np)';
      X = reshape(X_work(g,:,:),ncat,np)';
      P0 = hyge2pdf(0,N,x,X);
      P1 = hyge2pdf(1,N,x,X);
    else
      mu = reshape(mu_work(g,:,:),ncat,np)';
      P0 = binopdf(0,N,mu);
      P1 = binopdf(1,N,mu);
    end
    
    if any(P0(:)<0 | P0(:)>1 | P1(:)<0 | P1(:)>1)
      fprintf('WARNING:  invalid probabilities at some entries\n');
      keyboard;
    end
    
    % determine each patient's priority order of categories (according to P1)
    % left column of "priority" = least extreme category of mutation
    % right column of "priority" = most extreme category of mutation
    [tmp priority] = sort(P1,2,'descend');
    % sort the P arrays to put the columns in least->most priority order
    shft = (priority - repmat([1:ncat],np,1));
    map = reshape(1:(np*ncat),np,ncat);
    newmap = map + shft*np;
    P0 = P0(newmap);
    P1 = P1(newmap);
    
    P2 = 1-(P0+P1);  % note, P2 should really be called "P2_or_more"
    P2(P2<0) = 0;

    % STEP 2
    % for each sample, compute probability that it would have been of each (2-dimensional) degree.
    % degree=(d1,d2), where d=0 (no mut) ..... ncat (most extreme mut)
    % d1 is the MOST extreme mutation (or no mutation)
    % d2 is the SECOND MOST extreme mutation (or no mutation)
    % d1 can be 0-ncat; d2 can be 0-d1
    
    Pdeg = zeros(np,ncat+1,ncat+1);
    for d1=0:ncat, for d2=0:d1
      % has to have 0 in any/all categories > d1
      p = prod(P0(:,d1+1:end),2);
      if (d1>0)  % and (if d1>0)
        if (d1==d2)
          % if d1==d2, has to have 2+ in category d1
          p = p .* P2(:,d1);
        else
          % else:   has to have exactly 1 in category d1
          %         has to be clear in any/all categories (d2+1) to (d1-1)
          %         and (if d2>0) have (1 or 2+) in category d2
          p = p .* P1(:,d1);
          p = p .* prod(P0(:,d2+1:d1-1),2);
          if (d2>0)
            p = p .* (P1(:,d2)+P2(:,d2));
          end
        end
      end
      Pdeg(:,d1+1,d2+1) = p;
    end,end
    % sum of Pdeg(p,:,:) is confirmed to be 1 +- eps  (eps=2.2204e-16)

    %% STEP 2a: calculate score for a sample being of each possible degree
    %% (uses new style, where score = -log10 probability of having a mutation in that category
    %% (zero score for no mutations)
    Sdeg = zeros(np,ncat+1,ncat+1);
    for d1=1:ncat, for d2=0:d1
      if d1==d2
        p = P2(:,d1);
      else
        if d2>0
          p = P1(:,d1).*P1(:,d2);
        else 
          p = P1(:,d1);
        end
      end
      Sdeg(:,d1+1,d2+1) = -log10(p);
    end,end

    % CAPPING
    cap = -log10(P.projection_bin_cap);
    Sdeg(Sdeg>cap) = cap;

    % NO LANDFILLING


    % NULL BOOST
    if P.null_boost_factor>1
      null_boost_score = log10(P.null_boost_factor);
      nullcat = max(priority(:));  % assume null is the LAST category
      priority2 = [zeros(np,1) priority];
      Sdeg(priority2==nullcat) = Sdeg(priority2==nullcat) + null_boost_score;
      % TO DO:  Add double boost for two nulls
    end

    % STEP 3
    % determine actual (two-dimensional) degree and score for each sample
    % sum scores to get score_obs for gene
    
    degree = zeros(np,2);
    score_obs = 0;
    for p = 1:np
      i = 1;
      for d = ncat:-1:1
        c = priority(p,d);
        if i==1
          if n(p,c)>=2
            degree(p,:) = [d d];
            i = 3;
          elseif n(p,c)==1
            degree(p,i) = d;
            i=i+1;
          end
        elseif i==2
          if n(p,c)>=1
            degree(p,i) = d;
            i=i+1;
          end
        else % i>2: done
          break
        end
      end
      score_sample = Sdeg(p,degree(p,1)+1,degree(p,2)+1);
      score_obs = score_obs + score_sample;
    end

    % impose "tolerance" cutoff by decreasing "observed" score
    if P.mutsig_projection_tolerance_cutoff>0
      score_obs = score_obs * (1-P.mutsig_projection_tolerance_cutoff);
    end

    if score_obs<0
      fprintf('What?? score_obs<0???\n');
      score(g) = 0;
      pval(g) = 1;
    elseif isnan(score_obs)
      fprintf('What?? score_obs is nan?\n');
      keyboard
      score(g) = 0;
      pval(g) = 1;
    elseif score_obs==0
      score(g) = 0;
      pval(g) = 1;
    elseif isinf(score_obs)
      pval(g) = 0;
      score(g) = inf;
    else
      % STEP 4
      % compute P value for gene by convolutions

      % set up histogram
      %    -- choose binsize to be an even division of score_obs
      %       so that all probability in the last bin is *less* than the observed score
      minbins = P.convolution_minbins;
      if score_obs > minbins
        binsize = 1;
      else % even for very low-scoring genes, still make sure we use at least minbins bins
        binsize = score_obs / minbins;
      end
      numbins = ceil(score_obs / binsize);
      binsize = score_obs / numbins;
      if isinf(numbins)||isnan(numbins), fprintf('what?\n'); keyboard; continue;end

      H = zeros(numbins,1);
      H(1) = 1;  % initial condition: all probability is in first bin

      % sequential convolution
      offset = min(numbins, round(Sdeg/binsize));

      ncols = (ncat+1)*(ncat+2)/2;
      newH = zeros(numbins,ncols);
      for p = 1:np
        newH(:) = 0;
        col=1;
        for d1=0:ncat, for d2=0:d1
            o = offset(p,d1+1,d2+1);
            newH(o+1:end,col) = Pdeg(p,d1+1,d2+1) .* H(1:end-o);
            col=col+1;
        end,end
        H = sum(newH,2);
      end
      
      pval(g) = 1-sum(H);
      if pval(g)<0, pval(g) = 0; end
      score(g) = score_obs;

      if ~silence & ~debug
        if debug || (pval(g)<1e-4 || length(gta)<100)
          fprintf('%d\t%s\t%d\n', g, gname, pval(g));
        end
      end

%      if debug,       keyboard,      end

    end % next gene
    
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%% end of projection_2D method
  

case 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  BESTCAT METHOD
%
%  (1) For each category, calculate a p-value based on convoluting the patients.
%  (2) Then, take the best category (lowest p-value)
%  (3) Impose Bonferroni correction
%  (4) Assign that as the gene's p-value
%
%  Note: To consider "total" as a category, it must be explicitly included
%        as a column in the input matrix.
%
%  implemented 08-25-2012
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Using Bestcat method (for each category, convolute patients; then take best category.)\n');

% First check to see if "Total" is included as a category
tmp1 = sum(n_work(:,1:end-1,:),2);
tmp2 = n_work(:,end,:);
if all(tmp1(:)==tmp2(:))
  fprintf('\tVerified that Total has been included as final column of input matrices.\n');
else
  fprintf('\tWARNING:   TOTAL HAS NOT BEEN INCLUDED AS FINAL COLUMN OF INPUT MATRICES.\n');
end
clear tmp1 tmp2

silence = false;
debug = false;

fprintf('gene ');

for g=gta

  if mod(g,1000)==0, fprintf('%d/%d ', g, ng); end

  gname = P.gene_names{g};

  % save time if gene has no mutations
  if fullsum(n_work(g,1:ncat,:))==0
    pval(g) = 1;
    continue
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % TREAT EACH CATEGORY SEPARATELY
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  p_cat = nan(ncat,1);

  for c=1:ncat

    % save time if category has no mutations
    if sum(n_work(g,c,:))==0
      p_cat(c) = 1;
      continue
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 1
    % calculate point probability of the observation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Pobs = 1;
    for p=1:np
      N  = N_work(g,c,p);
      n  = n_work(g,c,p);
      
      if using_betabinom
        x = x_work(g,c,p);
        X = X_work(g,c,p);
        Pi = hyge2pdf(n,N,x,X);
      else
        mu = mu_work(g,c,p);
        Pi = exp(ln_nchoosek(N,n))*(mu^n)*((1-mu)^(N-n));
        if isinf(Pi)||isnan(Pi), Pi = binopdf(n,N,mu); end
      end
      Pobs = Pobs * Pi;
    end
    
    % impose minimum value
    cutoff = 1e-100;
    if Pobs<cutoff, Pobs = cutoff; end
    if Pobs>1, error('Pobs>1'); end
    
    if Pobs>0
      score_obs = -log10(Pobs);
    else
      %%% (should no longer happen, because we're using a cutoff of 1e-100
      keyboard
      fprintf('Pobs of zero in Bestcat for gene %s category %d!', gname,c);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 2
    % calculate P value of the observation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 2a
    % prepare table of piecewise probabilities
    % by patient and number of mutations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    MAXMUTS = P.convolution_maxmuts;
    Ppiece = zeros(np,MAXMUTS+1);    % second index is n+1 (to accommodate n=0)
    nmax = zeros(np,1);
    abort_flag = false;
    for p=1:np
      if abort_flag, break; end
      N = N_work(g,c,p);
      if N==0, continue; end   % if zero coverage, do not consider making any mutations
      if using_betabinom
        x = x_work(g,c,p);
        X = X_work(g,c,p);
      else
        mu = mu_work(g,c,p);
      end
      for n=0:MAXMUTS+1
        if abort_flag, break; end
        if using_betabinom
          pi = hyge2pdf(n,N,x,X);
        else
          pi = poisspdf(n,N*mu);
        end
        Ppiece(p,n+1) = pi;
        if pi < Pobs    % this individual score would alone "break the bank"
          n = n - 1;
          break
        elseif n > MAXMUTS
          fprintf('Exceeded MAXMUTS=%d with gene %d/%d (%s) category %d/%d pat %d/%d (%s): setting P=0\n',...
                  MAXMUTS,g,ng,gname,c,ncat,p,np,P.patient_names{p});
          p_cat(c) = 0;
          abort_flag = true;
          continue;
        end
      end
      nmax(p) = n;
    end
    maxnmax = max(nmax(:));
  
    if abort_flag, continue; end
    
    % compute scores and do landfilling
    Spiece = -log10(Ppiece);
    for p=1:np
      for n=0:nmax(p)-1
        if Spiece(p,n+1) >= Spiece(p,n+2)
            Spiece(p,n+1) = Spiece(p,c+2) - eps;
    end,end,end
    Spiece(:,1) = 0;  % no score for zero mutations

    % set up histogram
    % first: adjust binsize to make sure there is an integral number of bins per score_obs
    minbins = P.convolution_minbins;
    if score_obs > minbins
      binsize = 1;
    else % even for very low-scoring genes, still make sure we use at least minbins bins
      binsize = score_obs / minbins;
    end
    numbins = ceil(score_obs / binsize);
    binsize = score_obs / numbins;
    if isinf(numbins)||isnan(numbins), fprintf('what?\n'); keyboard; continue;end

    H = zeros(numbins,1);
    H(1) = 1;    % initial condition: all probability is in first bin (score=0, P=1)

    % sequential convolution
    offset = min(numbins,round(Spiece/binsize));
    newH = zeros(numbins,maxnmax+1);
    for p=1:np
      newH(:) = 0;
      for n=0:nmax(p)
        o = offset(p,n+1);
        newH(o+1:end,n+1) = Ppiece(p,n+1) .* H(1:end-o);
      end
      H = sum(newH(:,1:nmax(p)+1),2);
    end
    Pbulk = sum(H);
    p_cat(c) = 1 - Pbulk;
    
  end  % next category

  pval(g) = min(1,min(p_cat)*ncat);  % Bonferroni-corrected minimum across categories
  score(g) = -log10(pval(g));

  if ~silence
    fprintf('%d\t%s\t%d\n', g, gname, pval(g));
  end

end  % next gene

%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of Bestcat method


end % switch(METHOD)




%% finally, fix problematic p-values

pval(pval<0) = 0;
ntot = sum(sum(n_work,3),2);
Ntot = sum(sum(N_work,3),2);
pval(pval>1 | ntot==0 | Ntot==0) = 1;

fprintf ('Done\n');

