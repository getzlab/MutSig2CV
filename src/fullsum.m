function s = fullsum(m)
if length(m)==1, s = m;
else s = fullsum(sum(m)); end
end
