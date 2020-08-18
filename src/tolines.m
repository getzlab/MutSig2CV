function a = tolines(a)
a = split(a,char(10));
a = a(~cellfun('isempty',a));
