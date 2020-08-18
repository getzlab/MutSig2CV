function S = colx(n)

S = repmat({'col'},size(n));
for i=1:length(n)
  S{i} = [S{i} num2str(n(i))];
end
