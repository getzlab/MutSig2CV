function x = add_and_convert_simple_fieldnames(x,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'force_recalculate_maf_simplename_fields',false);
P = impose_default_value(P,'quiet',1);

if ~isfield(x,'start') && isfield(x,'pos'), x.start=x.pos; end
if ~isfield(x,'end') && isfield(x,'pos'), x.end=x.pos; end

try
  x = add_simple_fieldnames(x,P);
catch me
  disp(me.message);
  fprintf('Warning: add_simple_fieldnames failed.  Maybe fieldnames are already simple.\n');
end

if isfield(x,'patient_name') && isfield(x,'patient'),
  % if patient looks like it's really patient_idx, then replace it with patient_name
  tmp = str2double(x.patient);
  if mean(isnan(tmp))<0.05
    x.patient = x.patient_name;
  end
end

if isfield(x,'gene_name') && isfield(x,'gene'),
  tmp = str2double(x.gene);
  if mean(isnan(tmp))<0.05
    x.gene = x.gene_name;
  end
end

if ~isfield(x,'newbase'), x.newbase = find_newbase(x); end
if ~isfield(x,'start') && isfield(x,'pos'), x = rename_field(x,'pos','start'); end

%mandatory_flds = {'patient','chr','start','end','type','gene','classification','ref_allele','newbase'};
mandatory_flds = {'patient','chr','start','end','gene','ref_allele','newbase'}; % classification and newbase can now be imputed by enforce_effect
require_fields(x,mandatory_flds);

x.chr = convert_chr(x.chr,P);
x = make_numeric(x,{'start','end'});

if isfield(x,'firehose_patient_id'), x.patient = x.firehose_patient_id; end   % MAKE SURE



