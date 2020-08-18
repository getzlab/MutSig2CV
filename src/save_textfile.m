function save_textfile(t,filename)
% save_textfile(t,filename)
%
% Mike Lawrence 2008-10-15
out = fopen(filename,'wt');
fwrite(out,t,'char');
fclose(out);
