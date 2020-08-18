function X = load_struct_specify_string_cols(fname,string_cols,num_header_lines,lowercase_fieldnames)
% load_struct_specify_string_cols(fname,string_cols,num_header_lines)
%
% Mike Lawrence 2009-07-01

if ~exist('fname','var'), error('Must specify fname'); end
if ~exist('string_cols','var'), string_cols = []; end
if ~exist('num_header_lines','var'), num_header_lines = 1; end
if ~exist('lowercase_fieldnames','var'), lowercase_fieldnames = false; end

if ~exist(fname,'file'), error('%s not found',fname); end

numcols = get_colcount(fname,num_header_lines);

is_string = false(1,numcols);
is_string(string_cols) = true;
format = [];
for i=1:numcols
  if is_string(i), format = [format '%s'];
  else format = [format '%f']; end
end

%X = load_struct(fname,format,num_header_lines,lowercase_fieldnames);
P=[]; P.lowercase_fieldnames = lowercase_fieldnames;
X = load_struct(fname,format,char(9),num_header_lines,P);

