function p_out = fisher_combined_p(p_in)
% p_out = fisher_combined_p(p_in)
%
% p_in: each column is a set of p-values, assumed to be independent of each other
%       each row is a hypothesis

X2 = -2 * sum(log(p_in),2);
p_out = 1 - chi2cdf(X2,2*size(p_in,2));

