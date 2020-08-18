function numcols = get_colcount(fname,num_header_lines)
% Mike Lawrence 2009-10-16

if ~exist('num_header_lines','var'), num_header_lines = 1; end
if ~exist(fname,'file'), error('%s not found',fname); end

f = fopen(fname);
for i=1:num_header_lines+1; l = fgetl(f); end
numcols = sum(l==char(9))+1;
fclose(f);
