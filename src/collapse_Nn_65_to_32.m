function out = collapse_Nn_65_to_32(in)
% rows:  65 categories as in /xchip/tcga_scratch/lawrence/db/context65
% columns:  [N ->A ->C ->G ->T]
% Mike Lawrence 2010-01-27

if size(in,1)~=65, error('input must have 65 rows'); end
if size(in,2)~=5, error('input must have 5 columns (N A C G T)'); end

% (discard bottom row, "N")
out = collapse_Nn_64_by_strand(in(1:64,:));


