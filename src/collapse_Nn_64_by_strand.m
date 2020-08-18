function out = collapse_Nn_64_by_strand(in)
% rows:  64 categories as in /xchip/tcga_scratch/lawrence/db/allcateg/categs64.txt
% columns:  [N ->A ->C ->G ->T]
% Mike Lawrence 2009-12-11

if size(in,1)~=64, error('input must have 64 rows'); end
if size(in,2)~=5, error('input must have 5 columns (N A C G T)'); end

out = in(1:32,:,:,:,:);
compbase('ACGT') = 'TGCA';
X = generate_categ_context65_names();
for i=1:32
  oldname = X.name{i};
  newname = [compbase(oldname(1)) ' in ' compbase(oldname(end)) '_' compbase(oldname(end-2))];
  j = find(strcmp(newname,X.name));
  out(i,1,:,:,:) = in(i,1,:,:,:) + in(j,1,:,:,:);
  out(i,2,:,:,:) = in(i,2,:,:,:) + in(j,5,:,:,:);
  out(i,3,:,:,:) = in(i,3,:,:,:) + in(j,4,:,:,:);
  out(i,4,:,:,:) = in(i,4,:,:,:) + in(j,3,:,:,:);
  out(i,5,:,:,:) = in(i,5,:,:,:) + in(j,2,:,:,:);
end
