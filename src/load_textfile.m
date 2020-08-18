function t = load_textfile(filename)
in = fopen(filename);
t = fread(in,'uint8=>char')';
fclose(in);
