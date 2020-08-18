function x = parse_in(x,parsefield,pattern,newfields,numidx)
% x = parse_in(x,parsefield,pattern,newfields,numidx)

if ~exist('numidx','var'), numidx=[]; end

tmp = parse(getfield(x,parsefield),pattern,newfields,numidx);
x = merge_structs({x,tmp});
