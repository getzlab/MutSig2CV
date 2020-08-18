function keep = grepv(pattern,strings,flag)
% grepv(pattern,strings,flag)
%
% inverse grep
%
% Mike Lawrence 2010-02-04

if ~exist('flag','var'), flag=0; end

toss = grep(pattern,strings,1);
keep = setdiff((1:length(strings))',toss);
if ~flag, keep = strings(keep); end
