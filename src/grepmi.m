function match = grepmi(pat,str)

r = grepi(pat,str,1);
match = as_column(ismember(1:length(str),r));
