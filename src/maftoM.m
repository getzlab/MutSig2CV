function M = maftoM(maf,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'coding_only',false);

if isstruct(maf) && isfield(maf,'mut'), return; end  % already an M struct

if ischar(maf),
  mafname = maf;
  if exist(maf,'dir')
    maf = load2(mafname);
  elseif exist(maf,'file')
    maf = load_struct(mafname);
  end
end

if ~isfield(maf,'ttype'), maf.ttype = repmat({'---'},slength(maf),1); end
if ~isfield(maf,'patient') && isfield(maf,'Tumor_Sample_Barcode'), maf.patient = maf.Tumor_Sample_Barcode; end
if ~isfield(maf,'patient') && isfield(maf,'pat_idx'), maf.patient = maf.pat_idx; end
demand_fields(maf,{'patient','ttype'});

if ~isfield(maf,'chr') && isfield(maf,'Chromosome'), maf.chr=maf.Chromosome; end
if ~isfield(maf,'pos') && isfield(maf,'start'), maf.pos=maf.start; end
if ~isfield(maf,'pos') && isfield(maf,'Start_position'), maf.pos=maf.Start_position; end
if ~isfield(maf,'pos') && isfield(maf,'Start_Position'), maf.pos=maf.Start_Position; end

if isfield(maf,'chr'), maf.chr = convert_chr(maf.chr); end
if isfield(maf,'pos'), maf = make_numeric(maf,'pos'); end

if P.coding_only
  if isfield(maf,'type'), fld='type'; elseif isfield(maf,'Variant_Classification'), fld = 'Variant_Classification'; else fld = ''; end
%  fld = decell(grep('^(type|Variant_Classification)$',fieldnames(maf)));
  if isempty(fld), error('missing type information: can''t exclude noncoding.'); end
  is_coding = grepmi('silent|synonymous|missense|nonsense|splice|ins|del',maf.(fld));
  fprintf('Removing %d/%d noncoding mutations\n',sum(~is_coding),slength(maf));
  maf = reorder_struct_exclude(maf,~is_coding);
end


M=[];
if exist('mafname','var'), M.maf = mafname; end
M.mut = maf;
[tmp ui M.mut.pat_idx] = unique(M.mut.patient);
flds = {'patient','ttype','ttype0','source','repo','dataset'};
M.pat = keep_fields_if_exist(reorder_struct(M.mut,ui),flds);
M.pat = rename_field(M.pat,'patient','name');

M.pat.nmut = histc(M.mut.pat_idx,1:slength(M.pat));

if isfield(maf,'classification'), fld='classification'; elseif isfield(maf,'Variant_Type'), fld = 'Variant_Type'; else fld = ''; end
%fld = decell(grep('^(classification|Variant_Type)$',fieldnames(M.mut)));
if ~isempty(fld)
  midx = strcmp('SNP',M.mut.(fld)); M.pat.nsnp = histc(M.mut.pat_idx(midx),1:slength(M.pat));
  midx = strcmp('DNP',M.mut.(fld)); M.pat.ndnp = histc(M.mut.pat_idx(midx),1:slength(M.pat));
  midx = strcmp('DEL',M.mut.(fld)); M.pat.ndel = histc(M.mut.pat_idx(midx),1:slength(M.pat));
  midx = strcmp('INS',M.mut.(fld)); M.pat.nins = histc(M.mut.pat_idx(midx),1:slength(M.pat));
end

[M.ttype.name ui M.pat.ttype_idx] = unique(M.pat.ttype);
M.ttype.npat = histc(M.pat.ttype_idx,1:slength(M.ttype));
for i=1:slength(M.ttype), M.ttype.nmut(i,1) = sum(M.pat.nmut(M.pat.ttype_idx==i)); end

