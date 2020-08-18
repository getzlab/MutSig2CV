function pos2 = liftover(chr1,pos1,build1,build2,forceflag)
% pos2 = liftover(chr1,pos1,build1,build2,forceflag)

if ~exist('forceflag','var') forceflag=false; end

if nargin==2
  build1=18; build2=19;
  fprintf('Assuming you want hg18->hg19\n');
elseif nargin~=4
  fprintf('Need chr1,pos1,build1,build2\n');
end

if ~isnumeric(chr1), chr1 = convert_chr(chr1); end
if ~isnumeric(pos1), pos1 = str2double(pos1); end

if length(chr1)==1 && length(pos1)>1
  chr1 = chr1*ones(length(pos1),1);
end

if length(chr1)~=length(pos1), error('length(chr)~=length(pos)'); end

if isempty(pos1), pos2=[]; return; end

b1 = interpret_build(build1);
b2 = interpret_build(build2);
if ~((b1==18 && b2==19) || (b1==19 && b2==18) || (b1==38 && b2==19) || (b1==19 && b2==38))
  error('Unsupported combination of builds\n');
end

if b1==18 && b2==19
  fwb = '/cga/tcga-gsc/home/lawrence/db/hg18/liftOverToHg19/all.fwb';
elseif b1==19 && b2==18
  fwb = '/cga/tcga-gsc/home/lawrence/db/hg19/liftOverToHg18/all.fwb';
elseif b1==38 && b2==19
  fwb = '/cga/tcga-gsc/home/njharlen/db/hg38/liftOverToHg19/all.fwb';
elseif b1==19 && b2==38
  fwb = '/cga/tcga-gsc/home/njharlen/db/hg19/liftOverToHg38/all.fwb';
else
  error('wha?');
end

% make sure jar is on classpath
javaclasspath('/xchip/cga/reference/mutsig_params/FixedWidthBinary.jar')

F = org.broadinstitute.cga.tools.seq.FixedWidthBinary(fwb);

pos2 = double(F.get(chr1,pos1));

idx = find(pos2<1); %% positions that failed to map are returned as 0

if ~isempty(idx)
  fprintf('Note: %d/%d positions failed to map exactly.\n',length(idx),length(pos1));
  if ~forceflag
    fprintf('Returning these as NaN.  Use forceflag to force best match\n');
    pos2(idx) = nan;
  else
    fprintf('Forcing best match...\n');
    for i=1:length(idx), fprintf('%d/%d ',i, length(idx));
      xchr = chr1(idx(i));
      xpos = pos1(idx(i));
      ranges = [10,30,100,300,1000,3000,10000,100000,300000,1e6,3e6,10e6,30e6];
      rgcounts = zeros(size(ranges));
      for ri=1:length(ranges)
        range=ranges(ri);
        x1 = max(0,xpos-range):(xpos+range);
        x2 = double(F.get(xchr,x1(1),x1(end)));
        k = find(x2>=1,1);
        if ~isempty(k)
          diff = x1(k)-xpos;
          new = x2(k)+diff;
          pos2(idx(i))=new;
          rgcounts(ri)=rgcounts(ri)+1;
        end
      end
    end, fprintf('\n');
    fprintf('Ranges that were necessary:\n');
    pr(ranges,rgcounts,1:length(ranges));
  end
end

F.close();

