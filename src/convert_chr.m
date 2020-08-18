function C2 = convert_chr(C1,P)
%
% convert_chr(chromosome_list)
%
% converts text chromosome identifiers to numbers 1-24 (X=23, Y=24)
%
% Mike Lawrence 2008-05-01

if ~exist('P','var'), P=[]; end
if ischar(P)
  tmp=P;
  P=[];
  P.build = tmp;
end
P = impose_default_value(P,'build','');
P = impose_default_value(P,'quiet',0);

% CONVERSION METHODS AVAILABLE:
% if "build" is not specified, asssume human.
% if "build" is specified, but ReferenceInfoObj has not been initialized, use heuristics defined here.
% if "build" is specified AND ReferenceInfoObj has been initialized, use it to do the conversion.

RIO_available = false;

if isempty(P.build)
  assumed_human = true;
  ct = 24;
else
  try
    ct = ReferenceInfoObj.getMaxNum(P.build);
    RIO_available = true;
  catch me
    assumed_human = false;
    ct = get_chrcount(P.build);
  end
end

% CONVERSION

if isnumeric(C1)
  % do nothing
  C2=C1;

else

  if ~iscell(C1), C1={C1}; end
  if size(C1,1)==1, transpose_flag=true; C1=C1'; else transpose_flag=false; end

  [ccs cci ccj] = unique(C1);
  ccs_orig = ccs;

  if ~RIO_available
    % legacy behavior
    ccs = regexprep(ccs, '^([Cc][Hh][Rr])', '');
    ccs = regexprep(ccs, '^(Mt|MT)$','0');
    ccs = regexprep(ccs, '^[Mm]$', '0');
    
    idx = grep('^[XxYy]$',ccs,1);
    if ~isempty(idx)
      if assumed_human & ~P.quiet
        fprintf('convert_chr: assuming human for chrX/chrY\n');
      end
      ccs = regexprep(ccs, '^[Xx]$', num2str(ct-1));
      ccs = regexprep(ccs, '^[Yy]$', num2str(ct));
    end
    ccs_num = str2double(ccs);

  else
    % use ReferenceInfoObj
    ccs_num = nan(length(ccs),1);
    for i=1:length(ccs)
      try
        num = ReferenceInfoObj.getNum(ccs{i},P.build);
        if ~strcmp(ReferenceInfoObj.getUse(num,P.build),'D')
          num = inf;
        end
        ccs_num(i) = num;
      catch me
        % not found
      end
    end
  end

  idx = find(isnan(ccs_num));
  if ~isempty(idx) & ~P.quiet
    fprintf('WARNING: unknown chromosome(s):\n');
    disp(ccs_orig(idx));
  end
  idx = find(isinf(ccs_num));
  if ~isempty(idx) & ~P.quiet
    fprintf('WARNING: unused chromosome(s):\n');
    disp(ccs_orig(idx));
    ccs_num(idx) = nan;   % report as NaN
  end

  C2 = ccs_num(ccj);
  if transpose_flag, C2=C2'; end

end

idx = find(C2<0 | C2>ct);
if ~isempty(idx) & ~P.quiet
  fprintf('WARNING: chromosome number(s) out of range for this genome\n');
  disp(unique(C2(idx)));
end



