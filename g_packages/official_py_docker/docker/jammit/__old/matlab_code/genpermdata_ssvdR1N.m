function y = genpermdata_ssvdR1N(datacell,nperms,mode)
% y = genpermdata_ssvdR1N(datacell,nperms,mode)
% Input: mode = 0 (monte carlo); = 1 (permutations)
K = length(datacell);

y = cell(nperms,K);

for ii = 1:nperms
    
    for k = 1:K
        t = rowresamp(datacell{k},mode);
        y{ii,k} = t;
    end

end