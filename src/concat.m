function c = concat(strings,separator)
%
% concat(strings,separator)
%
% joins strings into one string, separated with the specified character/string
%
% Mike Lawrence 2008-05-01
% modified 2016-08-19 to work rowwise

if numel(strings)==length(strings)   % original behavior

  c='';
  for i=1:length(strings)
    if ischar(strings(i))
      c=[c strings(i)];
    elseif isnumeric(strings(i))
      c=[c num2str(strings(i))];
    elseif isnumeric(strings{i})
      c=[c num2str(strings{i})];
    elseif iscell(strings(i))
      c=[c strings{i}];
    else
      error('concat does not know how to handle that kind of input');
    end
    if i<length(strings)
      c=[c separator];
    end
  end

elseif length(size(strings))==2 && iscellstr(strings)     % modified to work rowwise

  nrow = size(strings,1); ncol=size(strings,2);
  c = cell(nrow,1);
  for j=1:nrow
    cc='';
    for i=1:ncol
      if ischar(strings(j,i))
        cc=[cc strings(j,i)];
      elseif isnumeric(strings(j,i))
        cc=[cc num2str(strings(j,i))];
      elseif isnumeric(strings{j,i})
        cc=[cc num2str(strings{j,i})];
      elseif iscell(strings(j,i))
        cc=[cc strings{j,i}];
      else
        error('concat does not know how to handle that kind of input');
      end
      if i<ncol
        cc=[cc separator];
      end
    end
    c{j}=cc;
  end

else
  error('concat does not know how to handle that kind of input');
end
