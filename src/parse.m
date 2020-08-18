function x = parse(S,r,f,numeric)
% x = parse(S,r,f)
%
% S = cell array of strings
% r = regexp to match
% f = fieldnames for output
% numeric = which columns to make numeric
%
% x = output struct with tokens from regexp stored in fields designated by f
%
% Mike Lawrence 2009-02-11 

if ~iscell(S), S = tolines(S); end
tokens = regexp(S,r,'tokens');

if ~iscell(f), f = {f}; end
x = tokens2struct(tokens,f);

if exist('numeric','var'), x = make_numeric(x,f(numeric)); end
