function varargout = parsein(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = parse_in(varargin{:});
else
  [varargout{1}] = parse_in(varargin{:});
end
