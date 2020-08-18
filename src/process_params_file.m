function P = process_params_file(P,fname)

if ischar(P)
  fname = P;
  P = [];
end

if ~isempty(fname)
  if ~exist(fname,'file')
    error('No such params file %s\n',fname);
  else
    tmp=[]; tmp.line = load_lines(fname);
    % try using whitespace delimiter
    tmp = parse_in(tmp,'line','^(\S+)\s*(\S*)$',{'key','value'});
    idx = find(strcmp(tmp.key,''));
    if ~isempty(idx)
      % try using tab delimiter
      tmp2 = parse_in(tmp,'line','^([^\t]+)\t*([^\t]*)$',{'key','value'});
      tmp.key(idx) = tmp2.key(idx);
      tmp.value(idx) = tmp2.value(idx);
    end
    for i=1:slength(tmp)
      key = tmp.key{i};
      value = tmp.value{i};
      if strcmpi(value,'true'), value=true; end
      if strcmpi(value,'false'), value=false; end
      if ischar(value)
        zzz = str2double(value);
        if ~isnan(zzz)
          value = zzz;
        else
          zzz = sscanf(value,'%d');
          if isnumeric(zzz) && ~isempty(zzz)
            value = zzz;
          end
        end
      end
      P = impose_default_value(P,key,value);
    end
  end
end
