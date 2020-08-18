function varargout = order_field_first(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = order_fields_first(varargin{:});
else
  [varargout{1}] = order_fields_first(varargin{:});
end
