function y = rowresamp(x,rflag)
% Usage: y = rowresamp(x,rflag)
% Input: x = pxn data matrix; rflag = 0 - bootstrap; 1 - permutation
% Output: y = pxn data matrix with resampled rows
% Description: This function resamples each row of x using permutations 
% or bootstrapping

if nargin < 2
    rflag = 0;
end

sz = size(x);
p = sz(1);
n = sz(2);
y = zeros(p,n);
%pop = x(:);
for ii = 1:p
    if rflag == 0
        y(ii,:) = randsample(x(ii,:),n,'true');  % Matlab function 
    else
        y(ii,:) = x(ii,randperm(n)); % Do permutation
    end
end
return;