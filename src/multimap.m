function [idx2,idx1] = multimap(S1,S2,flds1,flds2)
% [idx2,idx1] = multimap(S1,S2,flds)
% or
% [idx2,idx1] = multimap(S1,S2,flds1,flds2)
%
% analog of listmap() for structs, where mapping can be on basis of matching multiple fields
% 
% given structs S1 and S2,
%    maps S1 to S2 such that all fields specified in flds match
%
% flds1 and flds2 can both be specified, in case where the structs name the fields differently
%
% Mike Lawrence 2010-10-27

if ~exist('flds1','var'), flds1 = fieldnames(S1); end

if ~exist('flds2','var'), flds2=flds1; end

if length(flds1)~=length(flds2), error('length(flds1)~=length(flds2)'); end
nf = length(flds1);

demand_fields(S1,flds1);
demand_fields(S2,flds2);

method = 2;

if method==2

   % NEW METHOD

   totlen = slength(S1)+slength(S2);
   uj = nan(totlen,nf);
   
   for i=1:nf
     s1 = getfield(S1,flds1{i});
     s2 = getfield(S2,flds2{i});
     if (isnumeric(s1)|islogical(s1)) && (isnumeric(s2)|islogical(s2))
       % ok
     elseif iscellstr(s1) && iscellstr(s2)
       % ok
     else
       if isnumeric(s1)|islogical(s1), s1 = num2cellstr(s1); end
       if isnumeric(s2)|islogical(s2), s2 = num2cellstr(s2); end
     end
     [tmp tmp uj(:,i)] = unique([s1;s2]);
   end

   [tmp tmp vj] = unique(uj,'rows');
   st1 = 1;
   en1 = slength(S1);
   st2 = en1+1;
   en2 = totlen;

   idx2 = listmap(vj(st1:en1),vj(st2:en2));
   if nargout>1
    idx1 = listmap(vj(st2:en2),vj(st1:en1));
   end

elseif method==1

  % OLD METHOD, with stringsplice

  f1 = {};
  f2 = {};

  for i=1:nf
    q = getfield(S1,flds1{i});
    if isnumeric(q), q = num2cellstr(q); end
    f1 = [f1 q];
    q = getfield(S2,flds2{i});
    if isnumeric(q), q = num2cellstr(q); end
    f2 = [f2 q];
  end

  qf1 = stringsplice(f1,1,'###');
  qf2 = stringsplice(f2,1,'###');

  idx2 = listmap(qf1,qf2);
  if nargout>1
    idx1 = listmap(qf2,qf1);
  end


end





