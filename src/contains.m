function result = contains(string,substring)

if ~iscell(string), string={string}; end

if ~ischar(substring)
  error('second argument should be a string');
end

result = false(length(string),1);
for s=1:length(string)
  if length(substring)<=length(string{s})
    idx = 1:length(string{s})-length(substring)+1;
    for i=1:length(substring)
      idx = idx(find(string{s}(idx+i-1)==substring(i)));
      if isempty(idx), break; end
    end
    if ~isempty(idx), result(s) = true; end
  end
end  
