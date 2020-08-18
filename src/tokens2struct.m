function S = tokens2struct(T,fieldnames)
%
% Given the 'tokens' output of a regexp match,
% converts to a structure with the specified fieldnames
%
% e.g.
%
% >> T = regexp({'item=apple qty=2';'item=orange qty=4'},'item=(.*) qty=(.*)','tokens');
% >> S = tokens2struct(T,{'fruit','count'});
% >> look(S)
%
% [1]
%     fruit: apple
%     count: 2
%
% [2]
%     fruit: orange
%     count: 4
%
% Mike Lawrence 2009-02-04
%
nt = length(T);
nf = length(fieldnames);
S = [];

for f=1:nf
  F = cell(nt,1);
  for t=1:nt
    if isempty(T{t}), F{t} = '';
    else F{t} = T{t}{1}{f}; end
  end
  S = setfield(S,fieldnames{f},F);
end
