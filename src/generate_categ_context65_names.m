function categ_list = generate_categ_context65_names()

categ_list.num = (1:65)';
x = {'A';'C';'G';'T'};
y = {}; for i=1:length(x), y = [y; regexprep(x,'^(.*)$',[x{i} '_$1'])]; end
z = {}; for i=1:length(x), z = [z; regexprep(y,'^(.*)$',[x{i} ' in $1'])]; end
categ_list.name = [z;'any N'];

