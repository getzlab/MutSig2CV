function S = load_struct(varargin)
% load_struct(filename, format, separator_character, header_lines, <other_parameters>, P)

% loads a tab-delimited table from the specified file
% into a structure with field names based on the column headings
%
% piggybacks on read_table()
% all parameters are passed to read_table(),
% except last parameter if it is a struct (i.e. P struct)
% P struct can have following members:
%   P.lowercase_fieldnames: if true or 1, fieldnames are converted to lowercase
%
% WAS: load_struct(filename, format, header_lines, lowercase_fieldnames)
%   format is [] by default.
%   header_lines is 1 by default.
%   if header_lines is set to 0, field names are col1, col2, etc.
%   if lowercase_fieldnames is 1, fieldnames are converted to lowercase
%
% Mike Lawrence 2008-04-29
% fixed 2010-12-10 to give all its parameters to read_table
%   (except last parameter if a P struct)

%OLD:
%if ~exist('header_lines', 'var') || isempty(header_lines), header_lines = 1; end
%if ~exist('lowercase_fieldnames', 'var'),lowercase_fieldnames = false; end
%if ~exist('format', 'var'),format = [];end
%if ~exist(filename,'file'), error('%s not found',filename); end

% check parameters

args = varargin;

%%%%% DEPRECATED LEGACY CASE:
% X = load_struct(infile,'%s%s%f%f%s%s%s%s%s',0,[lowercase_fieldnames]);
if (length(args)==3 || length(args)==4) && ischar(args{2}) && isnumeric(args{3}) && ...
      (length(args)==3 || isnumeric(args{4}) || islogical(args{4}))
  fprintf('header_lines is now the fourth parameter... please update call to load_struct\n');
  args = [args(1:2) char(9) args(3:end)];
end

if length(args)<1
  error('filename required');
end
filename = args{1};
if ~ischar(filename)
  error('first parameter should be filename (character string)');
end
demand_file(filename);

% see if last parameter is a P struct
if isstruct(args{end})
  P = args{end};
  args = args(1:end-1);
else
  P = [];
end

P = impose_default_value(P,'lowercase_fieldnames',false);
P = impose_default_value(P,'ignore_poundsign_lines',true);
P = impose_default_value(P,'ignore_atsign_lines',true);
P = impose_default_value(P,'ignore_exclamationpoint_lines',false);

if length(args)>=2 
%  if ~ischar(args{2}) || length(args{2})<2 || args{2}(1)~='%'
%    error('second parameter should be format string, e.g. ''%s%f%f%f''');
%  end
end
if length(args)>=3
  if isempty(args{3})
    error('third parameter should be separator character, e.g. char(9)');
  end
  if isnumeric(args{3})
    error('header_lines is now the fourth parameter... please update call to load_struct');
  end
%  if ~ischar(args{3}) || length(args{3})~=1
%    error('third parameter should be separator character, e.g. char(9)');
%  end
end
if length(args)>=4
  if islogical(args{4})
    error('lowercase_fieldnames has been moved to P struct... please update call to load_struct');
  end
  if ~isnumeric(args{4})
    error('fourth parameter should be number of header lines');
  end
end

%OLD:
%table = read_table(filename, format, char(9), header_lines, 'whitespace', '\b\r', 'bufSize', 50000);

% HANDLE COMMENT LINES:
default_header_lines = 1;
%% see if table has comment lines at the beginning (start with #): if so, increment header_lines to skip them
n_comment_lines = 0;

if P.ignore_poundsign_lines
  f = fopen(args{1});
  while(true)
    x = fgetl(f);	
    if (isempty(x))||(x(1)=='#') % skip empty lines or lines starting with #
      n_comment_lines = n_comment_lines + 1;
      continue
    elseif strncmp(x,'Oncotator v',11)
      fprintf('Un-poundsigned Oncotator header detected and skipped.\n');
      n_comment_lines = n_comment_lines + 1;
      continue
    else
      break
    end
  end
  fclose(f);
end

%default_header_lines = 1 + n_comment_lines;

% HANDLE @ LINES (if we decide to remove them, useful for picard interval list headers): 
%% see if table has @ lines at the beginning (start with #): if so, increment header_lines to skip them

if P.ignore_atsign_lines
  f = fopen(args{1});
  while(true)
    x = fgetl(f);
    if (isempty(x))||(x(1)=='@') % skip empty lines or lines starting with #
      n_comment_lines = n_comment_lines + 1;
      continue
    elseif strncmp(x,'Oncotator v',11)
      fprintf('Un-poundsigned Oncotator header detected and skipped.\n');
      n_comment_lines = n_comment_lines + 1;
      continue
    else
      break
    end
  end
  fclose(f);
end


% HANDLE ! LINES (useful for reading Gene Ontology files)
if P.ignore_exclamationpoint_lines
  f = fopen(args{1});
  while(true)
    x = fgetl(f);
    if (isempty(x))||(x(1)=='!') % skip empty lines or lines starting with #
      n_comment_lines = n_comment_lines + 1;
      continue
    elseif strncmp(x,'Oncotator v',11)
      fprintf('Un-poundsigned Oncotator header detected and skipped.\n');
      n_comment_lines = n_comment_lines + 1;
      continue
    else
      break
    end
  end
  fclose(f);
end

default_header_lines = 1 + n_comment_lines;


%% default args
if length(args)==1, args = [args {''}]; end         % format string
if length(args)==2, args = [args {char(9)}]; end       % separator character
if length(args)==3, args = [args {default_header_lines}]; end          % number of header lines

%% ADDED BY PETAR -- we should be able to load a struct with no header, and still eliminate comment lines 

noheader_and_comments = 0;
if args{4} == 0 && default_header_lines > 1
  args{4} = default_header_lines-1; 
  noheader_and_comments = 1;
end 

% default whitespace
has_already = false;
for i=4:length(args)
  if ischar(args{i}) & strcmpi(args{i},'whitespace'), has_already=true; break; end
end
if ~has_already, args = [args {'whitespace'} {'\b\r'}]; end

% default bufSize - only usable pre-matlab 2014b
if verLessThan('matlab','8.4')
   has_already = false;
   for i=4:length(args)
      if ischar(args{i}) & strcmpi(args{i},'bufSize'), has_already=true; break; end
   end
   if ~has_already, args = [args {'bufSize'} {50000}]; end
end

% LOAD TABLE
try
  table = read_table(args{:});
  nf = length(table.dat);
catch me
  q = load_lines(args{1});
  if isempty(q)
    fprintf('\n%s is a blank file\n',args{1});
    table = [];
    table.dlm  = args{3};
    table.headers = {{}};
    table.dat = {{}};
    nf = 0;
  else
    disp(me);
    disp(me.message);
    error('Error loading struct file');
  end
end

if noheader_and_comments 
  table.headers = {}; 
end 
 
if isempty(table.headers)
  table.headers{1} = cell(nf,1);
  for f=1:nf
    table.headers{1}{f} = sprintf('col%d', f);
  end
end

% process header line
fields = table.headers{end};
if length(fields)~=nf
  fprintf('Header line has %d column headers instead of the expected %d:\n',length(fields),nf);
  fields{:}
  error('Unable to parse table header line.');
end

% remove illegal characters from column headings
% and convert to list of unique field names

if P.lowercase_fieldnames, fields = lower(fields); end
fields = regexprep(fields, '\W','');   % remove any characters except A-Z, a-z, 0-9, underscore
fields_orig = fields;
fields = genvarname(fields_orig);

% preserve "end", because it's only going to be a field name, not a variable name
for f=1:nf
  if strcmp(fields_orig{f}, 'end')
     fields{f} = 'end';
     break
  end
  if strcmp(fields_orig{f}, 'End')
     if P.lowercase_fieldnames
       fields{f} = 'end';
     else
       fields{f} = 'End';
     end
     break
  end
end

% convert table to structure

S = struct();
for f=1:nf
  S = setfield(S, fields{f}, table.dat{f});
end

% Added by Petar -- we should be able to select which fields we want to keep. 
% might be useful for loading huge maf files with "load_structs"

if isfield(P, 'fields_to_keep')

  if ~iscell(P.fields_to_keep) 
    error('"fields_to_keep" must be a cell array of strings!') 
  end 
  
  if sum(isnan(listmap(P.fields_to_keep, fieldnames(S)))) > 0
    error('"fields_to_keep" must be fields in the header of the struct!'); 
  end 

  S = keep_fields(S, P.fields_to_keep);

end 
