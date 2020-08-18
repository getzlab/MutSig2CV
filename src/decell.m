function x = decell(x)
if ~iscell(x) || length(x)~=1, error('improper use of decell(): expected a cell of length 1'); end
x = x{1};
