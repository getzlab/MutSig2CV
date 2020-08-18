function score = new_clustering_statistic(mutpos,genelength,metric,verbose_flag)
% score = new_clustering_statistic(mutpos,genelength)
%
% given a set of mutation positions (from 1 to genelength),
% calculates clustering statistic "score" (from 0 to 1)

if metric==1
  error('don''t call this function for metric1');
elseif metric==2
  d = diff(sort(mutpos));
  points = sum(1*(d<=20) + 1*(d<=12) + 2*(d<=6) + 4*(d<=2));
  maxpoints = 8 * (genelength-1);
  frac = points / maxpoints;
  score = frac .^ 0.5;  % to improve resolution at low end of range
elseif metric==2.5
  d = diff(sort(mutpos));
  points = sum(1*(d<=20) + 1*(d<=12) + 2*(d<=6) + 1*(d<=2) + 1*(d<=1) + 2*(d==0));
  maxpoints = 8 * (genelength-1);
  frac = points / maxpoints;
  score = frac .^ 0.5;  % to improve resolution at low end of range
elseif metric==2.6
  % Gaussian version

  if 0   % fit to a Gaussian
    d = 0:25;
    points25 = 1*(d<=20) + 1*(d<=12) + 2*(d<=6) + 1*(d<=2) + 1*(d<=1) + 2*(d==0);
    points2 = 1*(d<=20) + 1*(d<=12) + 2*(d<=6) + 4*(d<=2);
    gaussian6 = 120*normpdf(d,0,6);
    gaussian5 = 100*normpdf(d,0,5);
    gaussian45 = 90*normpdf(d,0,4.5);
    gaussian3 = 60*normpdf(d,0,3);
    plot(d,points2,'.-',d,points25,'.-',d,gaussian5,'.-','markersize',20)
    legend({'metric2','metric25','gaussian(5)'},'fontsize',20);
    xlabel('diff (bp)','fontsize',20); ylabel('score','fontsize',20);
    set(gcf,'color',[1 1 1]);
  end

  % gaussian1 = gaussian(6)
%  gaussian = [7.9788    7.8688    7.5477    7.0413    6.3890    5.6382    4.8394    4.0400    3.2802    2.5904    1.9895    1.4862    1.0798    0.7631 ...
%              0.5244    0.3506    0.2279    0.1441    0.0886    0.0530    0.0308    0.0175    0.0096    0.0051    0.0027    0.0014];

  % gaussian(5)
  gaussian = [7.9788    7.8209    7.3654    6.6645    5.7938    4.8394    3.8837    2.9945    2.2184    1.5790    1.0798    0.7095    0.4479    0.2717    0.1583...
              0.0886    0.0477    0.0246    0.0122    0.0058    0.0027    0.0012    0.0005    0.0002    0.0001    0.0000];

  % gaussian1.5 = gaussian(4.5)
%  gaussian = [7.9788    7.7842    7.2285    6.3890    5.3749    4.3038    3.2802    2.3796    1.6430    1.0798    0.6755    0.4022    0.2279    0.1229    0.0631];

  % gaussian2 = gaussian (3)
%  gaussian = [7.9788    7.5477    6.3890    4.8394    3.2802    1.9895    1.0798    0.5244    0.2279    0.0886    0.0308    0.0096    0.0027    0.0007];

  d = diff(sort(mutpos));
  d(d+1>length(gaussian))=[];
  points = sum(gaussian(d+1));
  maxpoints = gaussian(1) * (genelength-1);
  frac = points / maxpoints;
  score = frac .^ 0.5;  % to improve resolution at low end of range

elseif metric==2.7
  % try 1/x

  if 0   % fit to 1/x
    d = 0:25;
    points25 = 1*(d<=20) + 1*(d<=12) + 2*(d<=6) + 1*(d<=2) + 1*(d<=1) + 2*(d==0);
    points2 = 1*(d<=20) + 1*(d<=12) + 2*(d<=6) + 4*(d<=2);
    gaussian1 = 120*normpdf(d,0,6);
    gaussian15 = 90*normpdf(d,0,4.5);
    gaussian2 = 60*normpdf(d,0,3);
    hyp = 25 * (1./(d+3))-0.5;
    plot(d,points2,'.-',d,points25,'.-',d,gaussian1,'.-',d,gaussian2,'.-',d,gaussian15,'.-',d,hyp,'.-','markersize',20)
    legend({'metric2','metric2.5','gaussian(6)','gaussian(3)','gaussian(4.5)','hyp(24,3)'},'fontsize',20);
    xlabel('diff (bp)','fontsize',20); ylabel('score','fontsize',20);
    set(gcf,'color',[1 1 1]);
  end




elseif metric==3
  d = diff(sort(mutpos));
  points = sum(1*(d<=20) + 1*(d<=15) + 2*(d<=10) + 4*(d<=6));
  maxpoints = 8 * (genelength-1);
  frac = points / maxpoints;
  score = frac .^ 0.5;  % to improve resolution at low end of range
elseif metric==4
  d = diff(sort(mutpos));
  points = sum(1*(d<=24) + 2*(d<=12) + 3*(d<=6) + 5*(d<=3));
  maxpoints = 10 * (genelength-1);
  frac = points / maxpoints;
  score = frac .^ 0.1;  % to improve resolution at low end of range
elseif metric==5
  d = diff(sort(mutpos));
  points = sum(1*(d<=24) + 1*(d<=12) + 1*(d<=6) + 1*(d<=3));
  maxpoints = 4 * (genelength-1);
  frac = points / maxpoints;
  score = frac ;%.^ 0.5;  % to improve resolution at low end of range

elseif metric==100 || metric==101 || metric==102

  m = sort(mutpos);
  nm = length(m);
  L = genelength;

  % compute log-likelihood under uniform model
  uniform=[];
  uniform.LL = nm*log(1/L);

  % compute log-likelihood under nonuniform model
  %    find best region with increased mutation rate
  best=[];
  best.LL = -inf;
  for n1=1:nm
    for n2=n1:nm
      K = m(n2)-m(n1)+1;
      if metric==102 && K>10, continue; end  % limit region size 
      nk = n2-n1+1;
      nl = nm-nk;
      alpha = nl/nm;
      LL = nk*log((alpha/L)+((1-alpha)/K));
      if (nl>0), LL=LL+nl*log(alpha/L); end
      if LL>best.LL
        best.nk = nk;
        best.alpha = alpha;
        best.kstart = m(n1);
        best.kend = m(n2);
        best.kfrac = K/L;
        best.LL = LL;
  end,end,end

  if metric==101  % apply BIC
    best.LL = best.LL - log(nm);
  end

  score = best.LL - uniform.LL;
  score = score/100;

  if exist('verbose_flag','var') && verbose_flag==true
    best
    uniform
    score
  end

elseif metric==103

  m = mutpos;
  nm = length(m);
  L = genelength;

  % compute log-likelihood under uniform model
  uniform=[];
  uniform.LL = nm*log(1/L);

  % compute log-likelihood under nonuniform model
  %    find best 10bp region with increased mutation rate
  best=[];
  best.LL = -inf;
  K = 10;
  for n1=1:nm
    k1 = m(n1);
    k2 = k1+K-1;
    nk = sum(m>=k1 & m<=k2);
    nl = nm-nk;
    alpha = nl/nm;
    LL = nk*log((alpha/L)+((1-alpha)/K));
    if (nl>0), LL=LL+nl*log(alpha/L); end
    if LL>best.LL
      best.alpha = alpha;
      best.start = k1;
      best.end = k2;
      best.len = K/L;
      best.nmut = nk;
      best.LL = LL;
  end,end

  score = best.LL - uniform.LL;
  score = score/100;

  if exist('verbose_flag','var') && verbose_flag==true
    best
    uniform
    score
  end

elseif metric==200
  m = round(mutpos/3);
  nm = length(m);
  h=histc(m,1:max(m));
  thresh = max(2,nm * 0.05);
  score = sum(h(h>=thresh))/nm;

elseif metric==204   % decrease threshold from 5% to 2%
  m = round(mutpos/3);
  nm = length(m);
  h=histc(m,1:max(m));
  thresh = max(2,nm * 0.02);
  score = sum(h(h>=thresh))/nm;

elseif metric==201   % decrease threshold from 5% to 1%
  m = round(mutpos/3);
  nm = length(m);
  h=histc(m,1:max(m));
  thresh = max(2,nm * 0.01);
  score = sum(h(h>=thresh))/nm;

elseif metric==202   % replace round(mutpos/3) with running 12bp window
  m = sort(mutpos);
  nm = length(m);
  winsz = 12;
  thresh = max(2,nm * 0.01);
  ishot = false(nm,1);
  for winst=1:genelength-winsz
    inwin = (m>=winst & m<=(winst+winsz-1));
    if sum(inwin)>=thresh, ishot(inwin)=true; end
  end
  score = mean(ishot);

elseif metric==203   % replace round(mutpos/3) with running 36bp window
  m = sort(mutpos);
  nm = length(m);
  winsz = 36;
  thresh = max(2,nm * 0.01);
  ishot = false(nm,1);
  for winst=1:genelength-winsz
    inwin = (m>=winst & m<=(winst+winsz-1));
    if sum(inwin)>=thresh, ishot(inwin)=true; end
  end
  score = mean(ishot);

elseif metric==205   % replace round(mutpos/3) with running 6bp window
  m = sort(mutpos);
  nm = length(m);
  winsz = 6;
  thresh = max(2,nm * 0.02);
  ishot = false(nm,1);
  for winst=1:genelength-winsz
    inwin = (m>=winst & m<=(winst+winsz-1));
    if sum(inwin)>=thresh, ishot(inwin)=true; end
  end
  score = mean(ishot);

else
  error('unknown metric');
end








return

% visual scrambling
nx = 25;
for x=1:nx 
  subplot(ceil(sqrt(nx)),ceil(sqrt(nx)),x);
  plot(mutpos/3,randperm(length(mutpos)),'.')
  xlim([1 genelength/3]);
end













