function X = generate_1885_categ_names()

context = generate_categ_context65_names();
effect = get_effect29_categories_list();

X = cell(1885,1);
i=1;

for ci=1:slength(context)
  for ei=1:slength(effect)
    X{i} = [context.name{ci} ':' effect.name{ei}];
    i=i+1;
end,end

        
