function out = collapse_1885_to_192x5(in,dim)
% out = collapse_1885_to_192x5(in,dim)
%
% in = 1885 rows (corresponding to db/hg19/c65e29b)
%      (or 1885 elements in dimension "dim")
%      any number of other dimensions
%
% output = 192 rows (corresponding to reference/mutsig_params/192categs.txt)
%          5 columns (ncd/syn/mis/non/spl)
%          the other dimensions shifted forward one
%
% if dim>1, will singly expand the specified dimension in the output

if nargin~=1 && nargin~=2, error('wrong number of inputs'); end
if ~exist('dim','var'), dim=1; end

insz = size(in);
if insz(dim)~=1885, error('input must have 1885 elements in dimension %d',dim); end

% input assumes the ordering in db/hg19/c65e29b
%C1885 = load_struct('/cga/tcga-gsc/home/lawrence/db/hg19/c65e29b/categs.txt');
C1885 = generate_1885_categ_set();
C1885 = parse_in(C1885,'name','^(.+):(.+)$',{'context','effect'});
C1885 = parse_in(C1885,'context','^(.) in (.)_(.)$',{'base','left','right'});
C1885 = parse_in(C1885,'effect','^(...)/(...)/(...)$',{'eff1','eff2','eff3'});
idx=grep('noncoding',C1885.effect,1);z=repmat({'ncd'},length(idx),1);C1885.eff1(idx)=z;C1885.eff2(idx)=z;C1885.eff3(idx)=z;
idx=grep('splice',C1885.effect,1);z=repmat({'spl'},length(idx),1);C1885.eff1(idx)=z;C1885.eff2(idx)=z;C1885.eff3(idx)=z;
z=repmat({'---'},slength(C1885),1); C1885.effA=z; C1885.effC=z; C1885.effG=z; C1885.effT=z;
idx=find(strcmp('A',C1885.base));C1885.effC(idx)=C1885.eff1(idx);C1885.effG(idx)=C1885.eff2(idx);C1885.effT(idx)=C1885.eff3(idx);
idx=find(strcmp('C',C1885.base));C1885.effA(idx)=C1885.eff1(idx);C1885.effG(idx)=C1885.eff2(idx);C1885.effT(idx)=C1885.eff3(idx);
idx=find(strcmp('G',C1885.base));C1885.effA(idx)=C1885.eff1(idx);C1885.effC(idx)=C1885.eff2(idx);C1885.effT(idx)=C1885.eff3(idx);
idx=find(strcmp('T',C1885.base));C1885.effA(idx)=C1885.eff1(idx);C1885.effC(idx)=C1885.eff2(idx);C1885.effG(idx)=C1885.eff3(idx);
C1885 = rmfield_if_exist(C1885,{'eff1','eff2','eff3','effect','context'});

% output will use the default ordering of the 192 categories
%C192 = load_struct('/xchip/cga/reference/mutsig_params/192categs.txt');
C192 = generate_192_categ_set();
% output will use the canonical ordering of 'effect-based degrees'
effects = {'ncd','syn','mis','non','spl'};

%% DO THE MAPPING
newinsz = [insz(1:dim-1) 1885 1 insz(dim+1:end)]; % insert a singleton dimension to facilitate summing
in = reshape(in,newinsz);
outsz = [insz(1:dim-1) 192 5 insz(dim+1:end)];
out = zeros(outsz);
bases = 'ACGT';
for inrow=1:1885
  left = C1885.left{inrow};
  right = C1885.right{inrow};
  base = C1885.base{inrow};
  outrows = find(strcmp(C192.base,base) & strcmp(C192.left,left) & strcmp(C192.right,right));
  for newbase_idx=1:4
    newbase = bases(newbase_idx);
    if strcmp(newbase,base), continue; end
    outrow = outrows(strcmp(C192.newbase(outrows),newbase));
    if length(outrow)~=1, continue; end
    effect = C1885.(['eff' newbase]){inrow};
    effect_idx = find(strcmp(effects,effect));
    if length(effect_idx)~=1, continue; end
    outcol = effect_idx;
    idxstr = repmat(':,',1,length(insz)); idxstr(end)=[];
    outidxstr = [idxstr(1:dim*2-2) 'outrow,outcol' idxstr(dim*2:end)];
    inidxstr = [idxstr(1:dim*2-2) 'inrow,1' idxstr(dim*2:end)];
    cmd = ['out(' outidxstr ') = out(' outidxstr ') + sum(in(' inidxstr '),' num2str(dim) ');'];
    eval(cmd);
end,end
    




