function match = grepm(pat,str,varargin)
% returns boolean

if nargin>2
  fprintf('grepm: input arguments beyond first two are ignored\n');
end

r = grep(pat,str,1);
match = ismember((1:length(str))',r);
