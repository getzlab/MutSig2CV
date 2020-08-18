function q = remove_ansi(q)
q = regexprep(q,[char(27) '\[[^m]*m'],'');
