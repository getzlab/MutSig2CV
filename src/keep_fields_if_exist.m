function varargout = keep_fields_if_exist(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = keep_fields_that_exist(varargin{:});
else
  [varargout{1}] = keep_fields_that_exist(varargin{:});
end
