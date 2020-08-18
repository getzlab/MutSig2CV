function names = generate_192_categ_names

bases = 'ACGT';
names = cell(192,1);
i=1;
for from=1:4
  for left=1:4
    for right=1:4
      for to=1:4
        if from==to, continue; end
        names{i,1} = [bases(left) '(' bases(from) '->' bases(to) ')' bases(right)];
        i=i+1;
end,end,end,end
