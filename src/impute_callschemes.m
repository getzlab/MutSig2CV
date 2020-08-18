function M = impute_callschemes(M,P)

if ~exist('P','var'), P=[]; end
if isfield(P,'callschemes_all_WGS'), error('please use P.force_all_wgs_callschemes'); end
P = impose_default_value(P,'force_all_wgs_callschemes',false);
if P.force_all_wgs_callschemes, fprintf('Will assume all WGS callschemes.\n'); end

if isfield(M,'pat') && ~isfield(M,'patient')
  M = rename_field(M,'pat','patient');
  rename_pat_field_flag = true;
end

if P.force_all_wgs_callschemes
  M.patient.callscheme = 3*ones(slength(M.patient),1);
else
  if ~isfield(M.mut,'start')&&isfield(M.mut,'pos'), M.mut.start=M.mut.pos; end
  if ~isfield(M.mut,'end')&&isfield(M.mut,'pos'), M.mut.end=M.mut.pos; end

  % make sure M.mut has the necessary fieldnames
  M.mut = add_and_convert_simple_fieldnames(M.mut);
  M.mut = add_helper_is_fields(M.mut);

  if ~isfield(M,'patient')
    [M.patient.name tmp M.mut.pat_idx] = unique(M.mut.patient);
    M.np = slength(M.patient);
  end

  fprintf('Imputing callschemes\n');
  M.patient.nmut = as_column(histc(M.mut.pat_idx,1:slength(M.patient)));
  M.patient.n_coding = as_column(histc(M.mut.pat_idx(M.mut.is_coding),1:slength(M.patient)));
  M.patient.n_flank = as_column(histc(M.mut.pat_idx(M.mut.is_flank),1:slength(M.patient)));
  M.patient.n_coding_nonsilent = as_column(histc(M.mut.pat_idx(M.mut.is_coding & ~M.mut.is_silent),1:slength(M.patient)));
  M.patient.fracflank = M.patient.n_flank ./ M.patient.nmut;
  % 0=coding only;  1=exome+100bp flanks;  2=all capture (no interval list);  3=all genome (WGS)

  % first, estimate each patient individually, leaving NaN where too few mutations to know.
  M.patient.callscheme = nan(slength(M.patient),1);
  M.patient.callscheme(M.patient.nmut>=50 & M.patient.fracflank<0.03) = 0;
  M.patient.callscheme(M.patient.n_flank>=5 & M.patient.fracflank>=0.03) = 1;
  M.patient.callscheme(M.patient.n_flank>=5 & M.patient.fracflank>=0.3) = 2;
  M.patient.callscheme(M.patient.n_coding>=30 & M.patient.n_flank>=30 & M.patient.fracflank>=0.9) = 3;
  
  % if mutations are too few, then assume it's coding only
  M.patient.callscheme(isnan(M.patient.callscheme)) = 0;
  
  % if most patients are WGS, then call the rest of the patients WGS
  frac_wgs = mean(M.patient.callscheme==3);
  if frac_wgs>=0.8
    def=3;
  else % otherwise estimate the default capture callscheme, by summing together all NON-WGS SAMPLES
    pidx = find(M.patient.callscheme~=3);
    nmut_tot = sum(M.patient.nmut(pidx));
    n_flank_tot = sum(M.patient.n_flank(pidx));
    n_coding_tot = sum(M.patient.n_coding(pidx));
    fracflank_tot = n_flank_tot / n_coding_tot;
    def=1;  % this is the a priori most likely scheme given the current state of the pipeline
    if nmut_tot>=50 && fracflank_tot<0.03
      def=0;
    end
    if n_flank_tot>=5 && fracflank_tot>=0.03
      def=1;
    end
    if n_flank_tot>=10 && fracflank_tot>=0.3
      def=2;
    end
    if n_coding_tot>=30 && n_flank_tot>=30 && fracflank_tot>=0.9
      def=3;
    end
    % impose the default callscheme on the samples that had few mutations to know
    M.patient.callscheme(isnan(M.patient.callscheme))=def;
  end
  
  % if most capture patients are the same scheme, then override contrary call for patients with <100 muts
  cappats = find(M.patient.callscheme<3);
  if ~isempty(cappats)
    [md nummode] = mode(M.patient.callscheme(cappats));
    fracmode = nummode/length(cappats);
    if fracmode>0.8 && fracmode<1.0
      switchover = cappats(M.patient.nmut(cappats)<100 & M.patient.callscheme(cappats)~=md);
      if ~isempty(switchover)
        M.patient.callscheme(switchover) = md;
      end,end,end
end

callscheme_names = {'coding only';'exome+100bp flanks';'all capture (no interval list)';'all genome (WGS)'};
M.patient.callscheme_name = nansub(callscheme_names,M.patient.callscheme+1);

if exist('rename_pat_field_flag','var')
  M = rename_field(M,'patient','pat');
end




