function L = load_lines(fname)
X = load_textfile(fname);
L = text_to_lines(X);
