function A = num2strcell(a);

A = cell(length(a),1);
for i=1:length(a)
  A{i} = num2str(a(i));
end
