function categ_list = get_effect29_categories_list

changes = {'syn','mis','non'};
nc = length(changes);

n = nc^3;  % 3 possible newbases

categ_list = [];
categ_list.num = (1:n)';
categ_list.name = cell(n,1);

idx=0;
for i=1:nc
  for j=1:nc
    for k=1:nc
      idx=idx+1;
      categ_list.name{idx} = [changes{i} '/' changes{j} '/' changes{k}];
end,end,end

idx=idx+1;
categ_list.num(idx)=idx;
categ_list.name{idx} = 'splice-site';

idx=idx+1;
categ_list.num(idx)=idx;
categ_list.name{idx} = 'noncoding';

