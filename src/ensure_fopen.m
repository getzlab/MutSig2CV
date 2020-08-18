function result = ensure_fopen(varargin)
result = fopen(varargin{:});
if (result==-1) error('ERROR WRITING TO %s',varargin{1}); end

