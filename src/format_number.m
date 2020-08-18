function ss = format_number(values,sigfigs,width)
%
% format_number(value,sigfigs,width)
%
% formats the number specified in "value" to fit in
% a character field of "width" characters, while
% displaying "sigfigs" significant figures.
%
% if decimal expansion will fit, e.g. 0.0023, then this is used;
% otherwise scientific notation is used.
%
% Mike Lawrence 2008-05-20
%
% 2011-05-17 vectorized

if ~exist('sigfigs','var'), sigfigs=1; end
if ~exist('width','var'), width=4; end

if sigfigs<1, error('sigfigs must be at least 1'); end

if ~isnumeric(values), error('values must be numeric'); end
if ndims(values)>2, error('N-D case not implemented'); end

ni = size(values,1); nj = size(values,2);

if ni>1e5 || nj>1e5
  fprintf('format_number: table too big, using simpler method instead\n');
  ss = cellstr(num2str(values));

else

  sigfigs = sigfigs+1;

  ss = cell(ni,nj);
  for i=1:ni, for j=1:nj
      value = values(i,j);
      
      if value>0, pdz = ceil(-log10(value))-1;    % post-decimal zeroes
      elseif value<0, pdz = ceil(-log10(-value))-1;
      else, pdz=0; % value==0
      end

      s1 = sprintf(['%.' num2str(round(sigfigs)+pdz-(abs(value)<1)) 'f'], value);
      if str2double(s1)==0 && value~=0
        s1 = sprintf(['%.' num2str(1+round(sigfigs)+pdz-(abs(value)<1)) 'f'], value);
      end

      s2 = sprintf(['%.' num2str(round(sigfigs)-2) 'd'], value);
      
      if length(s1) <= width, s = s1;
      elseif length(s2) <= length(s1), s = s2;
      else, s = s1;
      end
      
      ss{i,j} = s;
  end,end

  if ni==1 && nj==1, ss = ss{1,1}; end   % to preserve original non-vectorized behavior

end
