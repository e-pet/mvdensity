function Y = soft_plus(X)
% SOFT-PLUS FUNCTION
%
% The soft-plus function y = log(1+exp(x)) is a smooth approximation to the
% rectifier y = max(0, x). 
% As exp(x) quickly exceeds the realm of standard floating point
% representations, we return y = x for x >= 40, for which the error is below 
% double precision: abs(log(1+exp(X))-X) < 1e-20. 
Y = zeros(size(X));
Y(X < 40) = log(1+exp(X(X < 40)));
Y(X >= 40) = X(X >= 40);
end