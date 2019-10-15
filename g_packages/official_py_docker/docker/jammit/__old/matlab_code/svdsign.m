function y = svdsign(data,k)
% Usage: y = svdsign(datamat,k)

sz = size(data);
nrows = sz(1);
ncols = sz(2);
%data = znrmmat(data,'rows','zscore');
[u,s,v] = svd(data,0);

kk = k-1;
datak = data;
for ii = 1:kk
   datak = datak-u(:,ii)*s(ii,ii)*v(:,ii)'; 
end
mudatak = mean(datak)';
cval = corr(mudatak,v(:,k));

if cval < 0
    y = -1;
else
    y = 1;
end
