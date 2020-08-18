function P = powerset(S)
  L = length(S);
  P = cell(2^L,1);
  for i=0:(2^L)-1
    idx = find(bitand(i,2.^(0:L-1)));
    P{i+1}=S(idx);
  end
end
