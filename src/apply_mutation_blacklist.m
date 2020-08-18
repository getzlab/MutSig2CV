function [x bad] = apply_mutation_blacklist(x,B,P)

if ischar(x), x=load_struct(x); end

if nargin==2 && isstruct(x) && isstruct(B)
  P=B;
  B=[];
end

if ~exist('x','var'), error('need at least mutation list'); end
bad = [];

if ~exist('P','var'), P=[]; end
if ~isfield(P,'build')
  fprintf('Assuming hg19\n');
%  fprintf('Is that OK? (dbcont/dbquit)\n'); keyboard
  P.build = 'hg19';
end

if ~exist('B','var') || isempty(B)
  def = '';
  if isfield(P,'build')
    b = interpret_build(P.build);
    if b==18
      def = '/xchip/cga/reference/mutsig_params/pancan_mutation_blacklist.v14.hg18.txt';
    elseif b==19
      def = '/xchip/cga/reference/mutsig_params/pancan_mutation_blacklist.v14.hg19.txt';
    end
  end
  if isempty(def)
    error('Don''t know what default blacklist to use');
  end
  fprintf('Loading default blacklist:\n\t%s\n',def);
  B = def;
end

if ischar(B)
  B = load_struct(B);
end

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'mutation_blacklist_match_fields','chr,start,newbase');

try
  B = add_simple_fieldnames(B);
catch me
end

try
  x = add_simple_fieldnames(x);
catch me
end

flds = split(P.mutation_blacklist_match_fields,',');

if ~isfield(x,'start') && isfield(x,'pos'), x.start=x.pos; end
if ~isfield(x,'end') && isfield(x,'pos'), x.end=x.pos; end

if ~all(isfield(x,flds)) || ~all(isfield(B,flds))
  fprintf('MAF is missing some fields specified in P.mutation_blacklist_match_fields, which are:\n');
  disp(flds);
else
  bi = multimap(x,B,flds);
  idx = find(~isnan(bi));
  if ~isempty(idx)
    fprintf('Omitting the following %d blacklisted mutations:\n', length(idx));
    pr(x,{'patient','gene','chr','start','newbase'},idx);
    if length(idx)>40, fprintf('Omitted %d blacklisted mutations\n', length(idx)); end
    %disp([x.patient(idx) x.gene(idx) num2cell([x.chr(idx) x.start(idx)]) x.newbase(idx)])
    if nargout>=2, bad = reorder_struct(x,idx); end
    x = reorder_struct_exclude(x,idx);
  end
end
