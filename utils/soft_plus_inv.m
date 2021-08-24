function Y = soft_plus_inv(X)
% INVERSE OF THE SOFT-PLUS FUNCTION
%
% The soft-plus function y = log(1+exp(x)) is a smooth approximation to the
% rectifier y = max(0, x). This function implements the inverse of that
% function, i.e., y = log(exp(x)-1).
% As exp(x) quickly exceeds the realm of standard floating point
% representations, we return y = x for x >= 40, for which the error is below 
% double precision: abs(log(exp(X)-1)-X) < 1e-20. 
Y = zeros(size(X));
Y(X < 40) = log(exp(X(X < 40))-1);
Y(X >= 40) = X(X >= 40);
end