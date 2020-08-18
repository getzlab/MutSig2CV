function R = split(string, delim)
% implementation of Perl "split"
% Mike Lawrence 2008-06-11

if ischar(string), string = {string}; charflag=true; else charflag=false; end

ns = length(string);

R = cell(ns,1);

for z=1:ns, if mod(z,10000)==0, fprintf('%d/%d ',z,ns); end

  dpos = [0 find(string{z}==delim) length(string{z})+1];
  nt = length(dpos)-1;
  R{z} = cell(nt,1);
  for t=1:nt
    R{z}{t} = string{z}(dpos(t)+1:dpos(t+1)-1);
  end

end, if ns>=10000, fprintf('\n'); end

if charflag, R = R{1}; end

return

% old implementation (treats multiple delimiters as one)

  t = regexp(string,['([^' delim ']*)'],'tokens');
  nt = size(t,2);
  R = cell(nt,1);
  for i=1:nt
    R{i} = t{i}{1};
  end
end
