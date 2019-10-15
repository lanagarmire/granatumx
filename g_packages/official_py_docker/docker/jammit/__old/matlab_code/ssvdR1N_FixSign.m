function y = ssvdR1N_FixSign(datacell,alpha) 
% Usage: y = ssvdR1N_FixSign(datacell,alpha)
%
% Input: datacell = cell array of K data matrices
%        alpha = l1 "sparsity" parameter   
% Output: y = K+1-dimensional cell array of K sparse loading 
%             vectors that model the common signal in y(K+1).
% Description: This version of the sparse SVD can handle an arbitrary 
%              matrices and fixes the sign problem inherent to the SVD.
% Get number of input data matrices of varargin
K = length(datacell);

% Get dimensions of the input data matrices
nrowvec = [];
for k = 1:K
    t = datacell{k};
    tt = size(t);
    tt = tt(1);
    nrowvec = [nrowvec tt];
end
nrowvec;

% ??? the tbegin and tend are not used
t = [1 nrowvec];
t0 = 0;
t1 = 0;
tbegin = [];
tend = [];
for k = 1:K 
    t0 = t0+t(k);
    tbegin = [tbegin t0];
    t1 = t1+t(k+1);
    tend = [tend t1];
end
rints = [tbegin' tend']';

% Form stack of input matrices
X = [];
for k = 1:K
    X = [X;datacell{k}];
end

% Compute initial solution based on the classical SVD
[t1, t2, t3] = svd(X,'econ');

% initial value of left singular vector
% ??? what is svdsign?
u0 = svdsign(X,1).*t1(:,1); %+ svdsign(X,2).*t1(:,2);

% ud is used to to decide convergence
ud = 1;
% number of iterations
iter_cnt = 0;

% Fix sign of u0
% mu = mean(X)';
% t31 = t3(:,1);
% corrmut31 = corr(mu,t31);
% if corrmut31 <= 0
%     u0 = -u0;
% end

% Compute sparse, rank-1 approximation of stacked matrix
while (ud > 0.0001) 
    iter_cnt = iter_cnt+1;    
    v =  X'*u0/sqrt(sum((X'*u0).^2)); %updating v 
    %u = sign(X*v).*(abs(X*v)-alpha);%updating u
    u = sign(X*v).*max(abs(X*v)-alpha,0);%updating u
    s = sqrt(sum(u.^2)); %singular value
    u = u/s;% normalizing u    
    ud = sqrt(sum((u0-u).^2));%ud=||u-u0||
    u0 = u;
    if iter_cnt > 5000 %5000 is the maximum number of iterations
        disp('Failed to converge! Increase the limit on the maximum number of iterations')
        break
    end    
end

% Define variable output cell
y = cell(K+1,1);

% Unwind stacked loading vector
for k = 1:K
    t = u([rints(1,k):rints(2,k)]);
    y{k} = t;
end

% Compute common dominant signal in row-space of stacked matrix X
V=sqrt(sum((X'*u).^2))*v;
y{K+1} = V;

return