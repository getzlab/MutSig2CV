function S = merge_structs(s1,s2)
%merge_structs takes two data structures s1 and s2 and returns a united
%structure S created by concating the fields of s2 to those of s1
%
%   S = merge_structs(s1,s2)
%   
%   s1 and s2 are non-empty data structures to unite.  merge_structs will
%   match the data in identically named fields of s1 and s2,
%   concatinating the values in s2 immediately after those in s1. 
%   If a field is missing in either s1 or s2, the corresponding values
%   are empty in S.
%
%   See also unite_Ds
%
%
%---
% $Id$
% $20080704 00:03$
% $cmermel$
% $1$

  
  n1 = size(s1,1);
  n2 = size(s2,1);
  S = struct();
  
  fields = union(fieldnames(s1),fieldnames(s2));
  
  str = 'struct(''';
  for i=1:length(fields)
    if isfield(s1,fields{i})
      xx1 = eval(['{s1.' fields{i} '}'])';
    else
      xx1 = cell(n1,1);
    end
    if isfield(s2,fields{i})
      xx2 = eval(['{s2.' fields{i} '}'])';
    else
      xx2 = cell(n2,1);
    end
    xx{:,i} = cat(1,xx1,xx2);
    str = cat(2,str,fields{i},''',','xx{',num2str(i),'}',',''');
  end
  
  str = [str(1:length(str)-2) ');'];
  eval(['S = ' str])
  
  
