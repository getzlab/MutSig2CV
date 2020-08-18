function x = add_simple_fieldnames(x,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'force_recalculate_maf_simplename_fields',false);

toflds = {'patient','type','classification','start','end',...
  'gene','chr','ref_allele','tum_allele1','tum_allele2','newbase'};
frflds = {'Tumor_Sample_Barcode','Variant_Classification','Variant_Type','Start_Position','End_Position',...
  'Hugo_Symbol','Chromosome','Reference_Allele','Tumor_seq_allele1','Tumor_seq_allele2','...get newbase'};  

if length(toflds)~=length(frflds), error('What?'); end

flds = fieldnames(x);

for i=1:length(toflds)
  to=[];
  if ~P.force_recalculate_maf_simplename_fields
    toidx = find(strcmpi(toflds{i},flds),1);
    if ~isempty(toidx)
      to = getfield(x,flds{toidx});
    end
  else
    % pretend there was nothing already named according to the "to" fieldname
    % (this flag is used in Firehose-deployed versions)
  end
  fr=[];
  if strcmpi(frflds{i},'...get newbase')
    fr = find_newbase(x);
  else
    fridx = find(strcmpi(frflds{i},flds),1);
    if ~isempty(fridx)
      fr = getfield(x,flds{fridx});
      if strcmpi(frflds{i},'Tumor_Sample_Barcode')
        fr = regexprep(fr,'^(.*)-Tumor$','$1');
      end
    end
  end
  if ~isempty(to) && ~isempty(fr)     % if it already has "to" AND "from", then fill in blanks in "to"
    if iscellstr(fr) && iscellstr(to)
      idx = find(strcmp('',to) & ~strcmp('',fr));
      to(idx) = fr(idx);
      x = setfield(x,flds{toidx},to);
    end
  elseif isempty(to) && ~isempty(fr)  % if it doesnt already have "to" and does have "from", then to=from
    x = setfield(x,toflds{i},fr);
  end
end   


