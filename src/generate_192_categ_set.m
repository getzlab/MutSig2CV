function C192 = generate_192_categ_set

C192=[]; C192.name = generate_192_categ_names();
C192 = parsein(C192,'name','^(.)\((.)\->(.)\)(.)$',{'left','base','newbase','right'});
C192.type = repmat({'point'},slength(C192),1);

