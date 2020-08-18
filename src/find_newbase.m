function nb = find_newbase(M)

%if isfield(M, 'newbase') & any(cellfun('isempty',M.newbase))
%M = reorder_struct(M, ~cellfun('isempty',M.newbase));
%end 
if isfield(M, 'newbase')
 reorder_struct(M, ~cellfun('isempty',M.newbase));
end 

%keyboard
if isfield(M,'newbase') && ~any(cellfun('isempty',M.newbase))
  nb = M.newbase;
else
%  keyboard
  if isfield(M,'TumorAllele')
    nb = M.TumorAllele;
  else
    if ~isfield(M,'ref_allele')
      M.ref_allele = M.Reference_Allele;
      M.tum_allele1 = M.Tumor_Seq_Allele1;
      M.tum_allele2 = M.Tumor_Seq_Allele2;
    end 
    %keyboard
    nb = M.tum_allele1;
    idx = find(strcmp(M.ref_allele,M.tum_allele1));
    nb(idx) = M.tum_allele2(idx);
  end
end
