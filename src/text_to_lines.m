function L = text_to_lines(X)

%X = regexprep(X,'\w+$','');

while ~isempty(X)
  if X(end)==10, X(end)=[];   % remove trailing blank lines
  else break; end
end

if isempty(X), L={};
else L = split(X,char(10)); end

