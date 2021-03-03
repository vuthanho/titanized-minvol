function [X] = simplexProj(Y,epsilon)
% this function projects each column x of a matrix X on the set defined by
%  {y in R^d s.t. y>=epsilon, sum_i y_i=1}
if nargin < 2
    epsilon = 0;
end
[n, ~] = size(Y);
b = ones(n,1);
lambda = max((cumsum(sort(Y.*b,1,'descend'))-1)./(cumsum(b.*b)));
X = max(epsilon, Y - lambda);
end