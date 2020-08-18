function s=reorder_struct_exclude(s,order)
if islogical(order), order = find(order); end
s = reorder_struct(s,setdiff(1:slength(s),order));
