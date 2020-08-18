function varargout = get_filedatenum(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = filedatenum(varargin{:});
else
  [varargout{1}] = filedatenum(varargin{:});
end
