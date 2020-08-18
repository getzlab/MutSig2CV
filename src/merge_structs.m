function out = merge_structs(varargin)
% merge_structs(s1,s2)
%
%    (two input parameters, each a struct)
%    calls merge_structs1.m  
%    written by Craig Mermel
%
% merge_stucts(slist)
%
%    (one input parameter, a cell array of structs)
%    calls merge_structs2.m
%    written by Mike Lawrence
%
%

if nargin==2
  out = merge_structs1(varargin{:});
elseif nargin==1 
  out = merge_structs2(varargin{:});
else
  error('Unknown input type for merge_structs');
end

