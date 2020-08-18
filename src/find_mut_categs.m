function [C,Sc,Sv] = find_mut_categs(Nn,P)
% Nn should be a 32x5 table as from collapse_Nn_64_by_strand
%   rows:  32 strand-collapsed categories
%   columns:  [N ->A ->C ->G ->T]
%
% Nn can also be an M struct with mut and cov fields (from load_all_mutation_data2.m)
%
% C is list of best category sets discovered
% Sc is stats of each best set of categories
% Sv is stats of "vanilla" category set of (CpG transit, other C:G transit, C:G transver, A:T mut)
%
% based on category_discovery.m, Mike Lawrence August 2008
% packaged into current function 2009-12-11
% modified to accept M struct, 2010-12-01

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'method','best');   %  "greedy" or "best"
P = impose_default_value(P,'max_k',5);
P = impose_default_value(P,'mutcategs_report_filename',[]);

if isstruct(Nn)     % input is an M-type struct
  M = Nn;
  m = reorder_struct(M.mut,strcmp('SNP',M.mut.classification));
  bases = {'A';'C';'G';'T'};
  m.to = listmap(m.newbase,bases);
  N = sum(M.cov.orig_cov)';
  if ~isfield(m,'context65'), m.context65 = m.context; end
  m = make_numeric(m,'context65');
  n = hist2d_fast(m.context65,m.to,1,65,1,4);
  Nn = collapse_Nn_65_to_32([N n]);
end

if size(Nn,1)~=32, error('input must have 32 rows'); end
if size(Nn,2)~=5, error('input must have 5 columns (N A C G T)'); end

orig_N = [repmat(Nn(1:16,1)',3,1);repmat(Nn(17:32,1)',3,1)];
orig_n = [Nn(1:16,[3 4 5])';Nn(17:32,[2 4 5])'];

if ~isempty(P.mutcategs_report_filename)
  report_fh = fopen(P.mutcategs_report_filename,'wt');
end

% reformat N and n (6x16) into N and n (4x4x2x3)

% dimensions:
%   (1)  5' base (1234=ACGT)
%   (2)  3' base (1234=ACGT)
%   (3)  "from" base (12=AC)
%   (4)  "to" outcome (1=transition, 2=flip_transversion, 3=skew_transversion)

n = zeros(4,4,2,3);
N = zeros(4,4,2,3);

for base5 = 1:4
  for base3 = 1:4
    for oldbase = 1:2
      for muttype = 1:3
        col = 4*(base5-1)+base3;
        switch oldbase
          case 1 % A
            switch muttype
              case 1, row = 2;  % A->G transition
              case 2, row = 3;  % A->T transversion(flip)
              case 3, row = 1;  % A->C transversion(skew)
            end
          case 2 % C
            switch muttype
              case 1, row = 6;  % C->T transition
              case 2, row = 5;  % C->G transversion(flip)
              case 3, row = 4;  % C->A transversion(skew)
            end
        end
        n(base5,base3,oldbase,muttype) = orig_n(row,col);
        N(base5,base3,oldbase,muttype) = orig_N(row,col);
end,end,end,end

Ntot = sum(N(:));
unsplit = {[1:4],[1:4],[1:2],[1:3]};

C = cell(P.max_k,1);   % output of category sets
Sc = cell(P.max_k,1);    % output of stats about category sets

% test "vanilla" set for comparison
if nargout>=3
  vanilla_leaf = [4687 4799 25343 29183];  % CpG transitions, other C:G transitions, C:G transversions, A:T mutations
  Sv = [];
  Sv.H = entropy_by_parts2(vanilla_leaf);
  subfprintf('VANILLA CATEGORY SET:\n');
  [tmp stats] = reportrule2(vanilla_leaf);
  Sv.std_tot = sum(stats.ci_high(:)-stats.rate)/1.98;
end

if ~strcmpi(P.method,'best'), fprintf('WARNING: methods other than "best" do not properly route output into C and H'); end

switch(lower(P.method))

  case 'greedy'
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % navigate breakdown tree (greedy algorithm)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  initial = {unsplit};
  subfprintf('k=1:\n\n');
  H_initial = entropy_by_parts(initial);
  stats = reportrule(initial);
  std_tot = sum(stats.ci_high(:)-stats.rate)/1.98;
  
  current = initial;
  for sz=2:P.max_k
    subfprintf('k=%d:  ', sz);
    H_current = entropy_by_parts(current);
    best_new = {};
    best_dH = 0;
    for t=1:length(current)
      rest = current(setdiff(1:length(current),t));
      parent = current{t};
      for d=1:4
        tobreak = parent{d};
        p = powerset(tobreak);
        minel = min(tobreak);
        for i=2:length(p)-1
          if ismember(minel,p{i})
            child1 = parent;
            child2 = parent;
            child1{d} = p{i};
            child2{d} = setdiff(tobreak, p{i});
            new_parts = [rest;{child1};{child2}];
            H_new = entropy_by_parts(new_parts);
            dH = H_new - H_current;
            if dH<best_dH
              best_dH = dH;
              best_new = new_parts;
    end,end,end,end,end
     
    stats = reportrule(best_new);
    std_tot = sum(stats.ci_high(:)-stats.rate)/1.98;
    current = best_new;
  end % next sz

    
  case 'best-old1'
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % find best possible category set of a given size
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  initial = {unsplit};
  leaves = {initial};

  for sz=1:P.max_k
    subfprintf('k=%d:  ', sz);
    
    if sz>1
      old_leaves = leaves;
      leaves = {};
      for l = 1:length(old_leaves)
        leaf = old_leaves{l};
        for t=1:length(leaf)
          rest = leaf(setdiff(1:length(leaf),t));
          parent = leaf{t};
          for d=1:4
            tobreak = parent{d};
            p = powerset(tobreak);
            minel = min(tobreak);
            for i=2:length(p)-1
              if ismember(minel,p{i})
                child1 = parent;
                child2 = parent;
                child1{d} = p{i};
                child2{d} = setdiff(tobreak, p{i});
                new_leaf = [rest;{child1};{child2}];
                leaves = [leaves;{new_leaf}];
      end,end,end,end,end
      leaves = remove_duplicate_leaves(leaves);
    end
    
    best_l = 0;
    best_H = Inf;
    for l=1:length(leaves)
      H = entropy_by_parts(leaves{l});
      if H<best_H
        best_l = l;
        best_H = H;
      end
    end
    
    % report results
    best_leaf = leaves{best_l};
    stats = reportrule(best_leaf);
    std_tot = sum(stats.ci_high(:)-stats.rate)/1.98;
  end % next sz


  case 'best-old2'
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % find best possible category set of a given size
  % OPTIMIZED VERSION
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
   initial = {unsplit};
   leaves = {initial};
   
   estimated_growth_factor = 30;   % actual split closer to 20, but let's err on generous side
   for sz=1:P.max_k
     subfprintf('k=%d:  ', sz);
     if sz>1
       old_leaves = leaves;
       typical_leaf = repmat({unsplit},sz,1);
       leaves = repmat(typical_leaf,length(old_leaves)*estimated_growth_factor,1);
       lastleaf = 0;
       for l = 1:length(old_leaves)
         leaf = old_leaves{l};
         for t=1:length(leaf)
           rest = leaf(setdiff(1:length(leaf),t));
           parent = leaf{t};
           for d=1:4
             tobreak = parent{d};
             p = powerset(tobreak);
             p = p(2:end-1);   % remove empty and full sets
             np = length(p);
             new_leaf = [rest;{parent};{parent}];
             new_leaves = repmat({new_leaf},np,1);
             child1i = length(new_leaf);
             child2i = child1i-1;
             for i=1:np
               new_leaves{i}{child1i}{d} = p{i};
               new_leaves{i}{child2i}{d} = setdiff(tobreak,p{i});
             end
             leaves(lastleaf+1:lastleaf+length(new_leaves)) = new_leaves;
             lastleaf = lastleaf + length(new_leaves);
           end,end,end
           leaves = remove_duplicate_leaves(leaves(1:lastleaf));
     end

     best_l = 0; best_H = Inf;
     for l=1:length(leaves)
       H = entropy_by_parts(leaves{l});
       if H<best_H, best_l = l; best_H = H; end
     end
     
     % report results
     best_leaf = leaves{best_l};
     stats = reportrule(best_leaf);
     std_tot = sum(stats.ci_high(:)-stats.rate)/1.98;
   end % next sz

   
  case 'best'
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % superfast version based on all-integer processing
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   unsplit = categ_to_int({[1 2 3 4] [1 2 3 4] [1 2] [1 2 3]});
   mask = uint16(15*16.^[0:3]);
   nybs = uint16((1:14)'*(16.^[0:3]));

   initial = unsplit;
   leaves = initial;
   for k=1:P.max_k
       subfprintf('k=%d\n',k);
       if k>1
         old_leaves = leaves;
         leaves = zeros(34^(k-1),k,'uint16');
         lastleaf = 0;
         for l=1:length(old_leaves)
           leaf = old_leaves(l,:);
           for c=1:k-1   % choose which category to split
             parent = leaf(c);
             for d=1:4   % choose which dimension to split along
               tobreak = bitand(parent,mask(d));
               tokeep = parent-tobreak;
               frags1 = bitand(tobreak,nybs(:,d));
               frags2 = tobreak-frags1;
               frags = [frags1 frags2];
               frags(~frags1|~frags2,:)=[];
               children = tokeep+frags;
               nc = size(children,1);
               leaves(lastleaf+1:lastleaf+nc,1:k-1) = repmat(leaf,nc,1);
               leaves(lastleaf+1:lastleaf+nc,[c k]) = children;
               lastleaf = lastleaf + nc;
       end,end,end
       leaves = unique(sort(leaves(1:lastleaf,:),2),'rows');
     end

     best_l = 0; best_H = Inf;

     if k<4
       for l=1:length(leaves)
         H = entropy_by_parts2(leaves(l,:));
         if H<best_H, best_l = l; best_H = H; end
       end
     else
       if k==4
         mx = (2^16)-1;
         lookup = nan(mx,1);
         for i=1:mx
           if any(bitget(i,[11 12 16])), continue; end  % those bits have no meaning
           c = int_to_categ(i);
           ns = sumpart(n,c);
           Ns = sumpart(N,c);
           f = Ns / Ntot;
           H_part = entropy(ns/Ns);
           lookup(i) = f*H_part;
         end
       end
       for l=1:length(leaves)
         H = 0;
         for i=1:k
           H = H + lookup(leaves(l,i));
         end  
         if H<best_H, best_l = l; best_H = H; end
       end
     end
           
     % report results
     best_leaf = leaves(best_l,:);
     [C{k} stats] = reportrule2(best_leaf);
     std_tot = sum(stats.ci_high(:)-stats.rate)/1.98;
     Sc{k}.H = best_H;
     Sc{k}.std_tot = std_tot;
   end
   
 otherwise
  error('Unknown P.method');
end

if ~isempty(P.mutcategs_report_filename)
  fclose(report_fh);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions

  function subfprintf(str,varargin)
    fprintf(str,varargin{:});
    if ~isempty(P.mutcategs_report_filename)
      fprintf(report_fh,str,varargin{:});
    end
  end

  function H = entropy_by_parts(parts)
    H = 0;
    for pp=1:length(parts)
      ns = sumpart(n,parts{pp});
      Ns = sumpart(N,parts{pp});
      f = Ns / Ntot;
      H_part = entropy(ns/Ns);
      H = H + f*H_part;
    end
  end

  function H = entropy_by_parts2(partsi)
    H = 0;
    for pp=1:length(partsi)
      part = int_to_categ(partsi(pp));
      ns = sumpart(n,part);
      Ns = sumpart(N,part);
      f = Ns / Ntot;
      H_part = entropy(ns/Ns);
      H = H + f*H_part;
    end
  end

  function x = sumpart(m,cut)
    x = m(cut{1},cut{2},cut{3},cut{4});
    x = sum(x(:));
  end

  function H = entropy(p)
    if p==0 || p==1
      p1 = 0;
      p2 = 0;
    elseif p<0 || p>1
      error('p must be between zero and one');
    else
      p1 = p*log2(p);
      p2 = (1-p)*log2(1-p);
    end
    H = -(p1+p2);
  end

  function s = convert_parts_to_rule(parts)
    bases='ACGT';
    change='tfs';
    np = length(parts);
    s = cell(np,1);
    for p=1:np
      part = parts{p};
      s{p} = [bases(part{1}) '[' bases(part{3}) '->' change(part{4}) ']' bases(part{2})];
    end
  end

  function s = rulestats(parts)
    np=length(parts);
    s.n = zeros(np,1);
    s.N = zeros(np,1);
    for p=1:np
      s.N(p) = sumpart(N,parts{p});
      s.n(p) = sumpart(n,parts{p});
    end
    [s.rate ci] = binofit(s.n,s.N);  
    s.ci_low = ci(:,1);
    s.ci_high = ci(:,2);
    rate_tot = fullsum(n)/fullsum(N);
    s.relrate = s.rate / rate_tot;
  end

  function stats = reportrule(parts)
    rules = convert_parts_to_rule(parts);
    rulenames = find_good_names_for_mutation_categories(rules);
    stats = rulestats(parts);
    [tmp ord] = sort(stats.relrate,'descend');
    for j=1:length(parts)
      i=ord(j);
      subfprintf('%25s   n %5.0f N %10.0f  rate %.2e (%sx)  ci %.2e to %.2e\n',...
              rulenames{i},stats.n(i),stats.N(i),stats.rate(i),...
              format_number(stats.relrate(i),3,4),stats.ci_low(i),stats.ci_high(i));
    end
    subfprintf('\n');
  end

  function [ck, stats] = reportrule2(partsi)
    for i=1:length(partsi), parts{i,1} = int_to_categ(partsi(i)); end
    rules = convert_parts_to_rule(parts);
    stats = rulestats(parts);
    tmp = parse(rules,'(\S+)\[(\S+)\->(\S+)\](\S+)',{'left','from','change','right'});
    ck = merge_structs({tmp,stats});
    ck.autoname = rules;
    ck.name = find_good_names_for_mutation_categories(ck.autoname);
    ck.type = repmat({'point'},slength(ck),1);
    [tmp ord] = sort(stats.relrate,'descend');
    for j=1:length(parts)
      i=ord(j);
      subfprintf('%25s   n %5.0f N %10.0f  rate %.2e (%sx)  ci %.2e to %.2e\n',...
              ck.name{i},stats.n(i),stats.N(i),stats.rate(i),...
              format_number(stats.relrate(i),3,4),stats.ci_low(i),stats.ci_high(i));
    end
    ck = reorder_struct(ck,ord);
    subfprintf('\n');
  end

  function LL = remove_duplicate_leaves(L)
    LL = L;
    % first convert all categories to integers
    nl = length(L);
    len = zeros(nl,1);
    for l=1:nl
      leaf = L{l};
      nc = length(leaf);
      ileaf = zeros(nc,1);
      for c=1:nc
        ileaf(c) = categ_to_int(L{l}{c});
      end
      L{l} = ileaf;
      len(l) = nc;
    end
    % compare leaves in groups of same-length leaves
    deleted = false(nl,1);
    [u ui uj] = unique(len);
    for i=1:length(u)
      idx = find(uj==i);
      x = zeros(length(idx),u(i));
      for j=1:length(idx)
        x(j,:) = L{idx(j)};
      end
      x = sort(x,2);   % sort each row
      [v vi vj] = unique(x, 'rows');
      % mark duplicates for deletion
      deleted(setdiff(idx,idx(vi))) = true;
    end
    % delete duplicates
    LL = LL(find(~deleted));
  end

  function i = categ_to_int(c)
    x = [c{1} c{2}+4 c{3}+8 c{4}+12]-1;
    i = sum(bitshift(1,x));
  end

  function c = int_to_categ(i)
    x = bitget(i,1:16);
    b = 1:4;
    c = { b(x(1:4)>0) b(x(5:8)>0) b(x(9:12)>0) b(x(13:16)>0) };
  end  

end % of main function

