function keep = grepi(pattern,strings,flag)
% grepi(pattern,strings,flag)
%
% case-insensitive grep
%
% Mike Lawrence 2010-02-12

if ~exist('flag','var'), flag=0; end

keep = grep(upper(pattern),upper(strings),1);
if ~flag, keep = strings(keep); end
