function k = map_categories_to_65(categs_txt)
% k = map_categories_to_65(categs_txt)
%
% categs_txt is path to file db/*/categs.txt listing the categories to be mapped
% collapses all categories containing "C in C_T" etc. to one category,
% with the 65th category being "any N".

if ischar(categs_txt)
  demand_file(categs_txt);
  C = load_struct(categs_txt);
else
  C = categs_txt;
end

demand_fields(C,{'num','name'});
if iscellstr(C.num)
  if any(strcmp(C.num,'0')), error('C has a category #0'); end
  if ~strcmp(C.num{1},'1'), error('C first category is not #1'); end
elseif isnumeric(C.num)
  if any(C.num==0), error('C has a category #0'); end
  if C.num(1)~=1, error('C first category is not #1'); end
else
  error('unknown format for "num"');
end

k = nan(slength(C),1);
base = 'ACGT';
cidx = 1;
for mid = 1:4
  for left = 1:4
    for right = 1:4
      str = [base(mid) ' in ' base(left) '_' base(right)];
      idx = grep(str,C.name,1);
      if isempty(idx), error('C lacks a "%s" category!',str); end
      if any(~isnan(k(idx)))
        disp(C.name(idx(~isnan(k(idx)))))
        error('The above categories would be double counted!');
      end
      k(idx) = cidx;
      cidx=cidx+1;
end,end,end

str = 'any N';
idx = grep(str,C.name,1);
if isempty(idx), error('C lacks a "%s" category!',str); end
k(idx) = cidx;

if any(isnan(k(idx)))
  disp(C.name(idx(isnan(k(idx)))));
  fprintf('WARNING: the above categories were not counted.\n');
end



