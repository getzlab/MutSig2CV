function S = make_numeric(S,varargin)
%
% for a given struct S,
%   performs str2double on all the fields specified
%
% Mike Lawrence 2009-02-24

fields = {};
for i=1:length(varargin)
  fields = [fields varargin{i}];
end

for i=1:length(fields)
  if isfield(S,fields{i})    
    x = getfield(S,fields{i});
    if strncmpi(fields{i},'chr',3)
      if any(strcmpi('chrx',x)|strcmpi('chry',x))
          fprintf('\n');
          fprintf('\t***********************************************\n');
          fprintf('\t*                                             *\n');
          fprintf('\t*                 W A R N I N G               *\n');
          fprintf('\t*                                             *\n');
          fprintf('\t*   You are removing chrX/Y mutations!        *\n');
          fprintf('\t*                                             *\n');
          fprintf('\t*   Please use convert_chr()                  *\n');
          fprintf('\t*          not make_numeric()                 *\n');
          fprintf('\t*                                             *\n');
          fprintf('\t***********************************************\n');
          fprintf('\n');
          fprintf('Pausing for 10 minutes.\n');
          pause(600);
      end
    end
    if isnumeric(x) || islogical(x)
      %    fprintf('Warning: %s is already numeric.\n',fields{i});
    else
      x = str2doubleq_wrapper(x);
      if mean(~isnan(x))<0.90, fprintf('Warning: field "%s" is <90%% numeric: returning some NaNs.\n',fields{i}); end
      S = setfield(S,fields{i},x);
    end
  else
    fprintf('No such field: %s\n', fields{i});
  end
end
