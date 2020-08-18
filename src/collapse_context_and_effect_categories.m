function Q = collapse_context_and_effect_categories(C,K)
% C is struct {num,name} from c65e29 track
% K is categ struct from category discovery
%
% returns Q,  rows = raw categories from C
%             columns = collapsed categories from K
%             pages = 4 newbases
%
%             each cell tells the sense of the mutation:
%                 0 = incompatible combination
%                 2 = synonymous
%                 3 = missense
%                 4 = nonsense / nonstop / splice

k = assign_65x4_to_categ_set(K);
c = map_categories_to_65(C);
w = nansub(k,c);

C = parse_in(C,'name','^([ACGT]).*:(.*)$',{'from','effect'});
C = parse_in(C,'effect','^(syn|mis|non)/(syn|mis|non)/(syn|mis|non)$',{'o1','o2','o3'});
s = [listmap(C.o1,{'syn','mis','non'}) listmap(C.o2,{'syn','mis','non'}) listmap(C.o3,{'syn','mis','non'})];
s(strcmp(C.effect,'splice-site'),:) = 3;

Q = zeros(size(w));
bases = 'ACGT';
for from=1:4
  rows = find(strcmp(C.from,bases(from)));
  scol = 1;
  for newbase=1:4
    if newbase==from, continue; end
    Q(rows,:,newbase) = bsxfun(@times,w(rows,:,newbase),s(rows,scol));
    scol=scol+1;
  end
end

%                 0 = incompatible combination
%                 1 = synonymous
%                 2 = missense
%                 3 = nonsense / nonstop / splice

Q(Q>0) = Q(Q>0) + 1;

%                 0 = incompatible combination
%                 2 = synonymous
%                 3 = missense
%                 4 = nonsense / nonstop / splice
