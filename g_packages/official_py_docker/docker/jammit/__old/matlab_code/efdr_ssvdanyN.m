function y = efdr_ssvdanyN(datacell,nsuppvec,alphavec,nperms,...,
    rownames,colnames)
%
% Usage: y = efdr_ssvdanyN(datacell,nsuppvec,alphavec,nperms,...,
%   rownames,colnames)
%
% efdr_ssvdR1N.m jointly analyzes an arbitrary number of input data matrices
% and produces low-dimensional signatures for each matrix composed of a 
% small number of variables (i.e., matrix rows) that in aggregate represents 
% a dominant signal in the row-space of the matrix stack as a 
% sparse linear model. Different solutions are computed on a user-defined 
% grid of "sparsity" parameters and a false discovery rate (FDR) is 
% computed for each parameter based on permutation testing. Specific 
% solutions are selected by the user based on FDR and signature size.
%
% Input variables:
%   datacell = cell array of K input data matrices.
%   nsuppvec = vector of K estimates of signal support for each matrix.
%   alphvec = monotonically increasing grid of L sparsity params
%   nperms = number of permutations for estimating FDR
%   rownames = cell array of row labels for each data matrix
%   colnames = cell array of column labels for each data matrix
%
% Current output structure: 
% y=struct('nrowvec',nrowvec,'fdrmat',fdrmat,'nsig0mat',nsig0mat,...
%     'fdrtab',fdrtab,'fdrtable',fdrtable,'fdrbestmat',fdrbestmat,...
%     'ybestcell',{ybestcell},'sigmatcell matrix.
%   alphavec = grid of sparsity parameters alpha for FDR table.
%   nperms = number of permutations for each alpha of FDR t',{sigmatcell},...
%     'sigindxcell',{sigindxcell},'siglabelscell',{siglabelscell},...
%     'sumkmat',sumkmat);
%
% Output variables:
%   nrowvec = vector of row dimensions for each input matrix
%   fdrmat = matrix of estimated FDRs as a function of alpha (rows)
%       for each data matrix (columns).
%   nsig0mat = number of significant variables found by JAMMIT for each
%       matrix as a function of a user-defined grid of alphas.    
%   fdrtab = FDR table as a matrix of real numbers w/ columns for 
%       model indexing and alphas
%   fdrtable = fdrtab in Matlab table format w/ variable names.
%   fdrbestmat = matrix of fdrs for type-specific and aggregate signatures
%       as a function of the models selected by the user. 
%   ybestcell = cell array of data for the models selected by the user.
%   sigmatcell = cell array of signature data matrices.
%   signindxcell = cell array of signature indices in each matrix
%   siglabelscell = cell arry of signature labels based on input data for
%       row and column labels.
%   sumkmat = matrix of sparse linear models of dominant signal
%       in the stacked matrix as a function of selected solutions (rows)
%       over N samples (columns). 
%
disp('Start JAMMIT analysis!!');
tstart = tic;

% Get number of input data matrices
K = length(datacell);
sz = size(datacell{1});
N = sz(2);

% Display key input variables
% datacell
% rownames
% colnames
% nsuppvec
% alphavec
% nperms

% Get row-dimension of the input data matrices
nrowvec = [];
for k = 1:K
    t = size(datacell{k});
    tt = t(1);
    nrowvec = [nrowvec tt];
end

% Compute pi0  based on nsuppvec for each matrix
pi0vec = [];
for k = 1:K
    t = 1 - (nsuppvec(k)/nrowvec(k));   
    pi0vec = [pi0vec t];
end

% Compute estimated pi0 for stacked matrix
pi0all = 1-((sum(nsuppvec))/(sum(nrowvec)));
pi0vec = [pi0vec pi0all];

% Get length of sparsity parameter vector
nalpha = length(alphavec);

% Initialize vector of # of significant variables
nsig0mat = -ones(nalpha,K+1);

% Initialize FDR table components
fdrmat = -ones(nalpha,K+1);


% Compute nperms permuted data sets 
mode = 1; %(1 = permutation; 0 = montecarlo)
permdata = genpermdata_ssvdR1N(datacell,nperms,mode);

% Loop for each alpha value
for ii = 1:nalpha
    ii      
    % Get alpha 
    alphaii = alphavec(ii);
        
    % Compute JAMMIT soluton for alphaii
    y = ssvdR1N_FixSign(datacell,alphaii);
    
    % Find number of signigicant genes in the input data set
    nsig0ii = -ones(1,K);
    for k = 1:K
        uk = y{k};
        nsig0ii(k) = length(find(abs(uk)>0));  
    end
    nsig0ii = [nsig0ii sum(nsig0ii)];
    nsig0mat(ii,:) = nsig0ii;
    
    % Compute l1 residuals based on selected alpha
    residcell = {};
    for k = 1:K
        t = datacell{k}-y{k}*y{K+1}';
        residcell = [residcell {t}];
    end
    
%     % Compute nperms permuted data sets 
%     mode = 1; %(1 = permutation; 0 = montecarlo)
%     permdata = genpermdata_ssvdR1N(datacell,nperms,mode); 
%     % permdata = genpermdata_ssvdR1N(residcell,nperms);
    
    % Begin looping through the nperms permuted data
    nsigpermii = -ones(nperms,K);
    for jj = 1:nperms
      
        % Get permuted data
        yjj = permdata(jj,:);
        
        % Compute sparse solution on permuted data for selected alpha
        ujj = ssvdR1N_FixSign(yjj,alphaii);
                
        % Count false positives
        for k = 1:K
            nsigpermii(jj,k) = length(find(abs(ujj{k})>0));
        end
        
    end
    nsigpermallii = sum(nsigpermii,2); %mean(sum(nsigpermii,2));
    nsigpermii = [nsigpermii nsigpermallii];
    
    % Compute the mean number of false positives over nperms
    musigpermii = -ones(1,K+1);
    for k = 1:K+1
        musigpermii(k) = mean(nsigpermii(:,k));
    end
    
    % Compute K+1 FDRs for sparsity parameter alphaii
    fdrmat(ii,:) = min(pi0vec.*(musigpermii./nsig0ii),1);  
end

% Add number of significant variables over all K data matrices
nrowvec = [nrowvec sum(nrowvec)];

% Initialize FDR table
if K > 1
    fdrtab = -ones(nalpha,2*(K+1));
else
    fdrtab = -ones(nalpha,2*K);
end
%fdrtab
% Fill FDR table
for ii = 1:nalpha 
    t0 = [];
    t1 = nsig0mat(ii,:);
    t2 = fdrmat(ii,:);
    
    if K > 1       
        for k = 1:K+1
            t0 = [t0 [t1(k) t2(k)]];
        end
    else
        t0 = [t0 [t1(k) t2(k)]];
    end
    %t0
    %length(t0)
    fdrtab(ii,:) = t0;    
end

% Augment FDR table with alpha values and row indices
fdrtab = [[1:nalpha]' fdrtab alphavec'];
%fdrtab

% Do FDR table
%fdrtab

% Build string of table variables for FDR table
str=[];
if K > 1
    kend = 2*K+2+1+1;
    for k = 1:kend
        if k < kend
            str = [str 'fdrtab(:,',num2str(k),'),']; 
        else
            str = [str 'fdrtab(:,',num2str(k),')']; 
        end
    end
else
    kend = 2*K+1+1;
    for k = 1:kend
        if k < kend
            str = [str 'fdrtab(:,',num2str(k),'),'];
        else
            str = [str 'fdrtab(:,',num2str(k),')'];
        end
    end
end
%str

% Build string of variable names for FDR table
strlabels=[];
if K>1
    kend = K+1+1+1;
    for k = 1:kend
        if k == 1
            strlabels=[strlabels '''solution'','];
        elseif k < kend-1
            nsigstr=['''nsig',num2str(k-1),''','];
            fdrstr=['''fdr',num2str(k-1),''','];
            strlabels=[strlabels [nsigstr fdrstr]];
        elseif k == kend-1
            nsigstr=['''nsigall'','];
            fdrstr=['''fdrall'','];
            strlabels=[strlabels [nsigstr fdrstr]];
        else k == kend
            strlabels=[strlabels '''alpha'''];
        end
    end
else
    kend = K+1+1;
    for k = 1:kend
       if k == 1
          strlabels=[strlabels '''solution'','];
       elseif k == 2
          nsigstr=['''nsig',num2str(k-1),''','];
          fdrstr=['''fdr',num2str(k-1),''','];
          strlabels=[strlabels [nsigstr fdrstr]]; 
       else
          strlabels=[strlabels '''alpha''']; 
       end
    end
end

disp('End FDR estimation!!!');
% str
% strlabels

% Build string to generate FDR table
strvar=['''VariableNames'',{',strlabels,'}'];
strtable=['table(',str,',',strvar,');'];
%strvar
%strtable

% Compute FDR table
strdotable=['fdrtable=',strtable];
%strdotable

eval(strdotable);

% Display FDR table
fdrtable

% Get user-selected solution

% Initialize variables for user-defined solution
pflag = 1;
count = 0;
fdrbestmat = [];
sumkmat = [];
ybestcell = {};
sigmatcell = {};
sigindxcell = {};
siglabelscell = {};

telapsed = toc(tstart)

% Enable user selection of JAMMIT solution based on FDR table
while pflag
    count = count +1;
    
    if count == 1
        prompt0 = ...
            'Select desired row # from FDR table (Enter -1 to exit): ';
        disp('');
        % Prompt for "optimal" line of the FDR table
        result = input(prompt0);
        
        if result < 0
            y = {};
            disp(...
                ['Exiting JAMMIT pipeline with ',num2str(count-1),' models selected!!']);
            return;
        end
    end
    
        % Get and store selected row of FDR table
    fdrbest = fdrtab(result,:)
    
    % Store selected FDR in fdrbestmat matrix
    fdrbestmat = [fdrbestmat;fdrbest];
    % Get optimal alpha value
    alphabest = fdrbest(end);
    % Compute JAMMIT solution based on user-selected alpha
    ybest = ssvdR1N_FixSign(datacell,alphabest);

    % Build cell array of best solutions
    ybestcell{length(ybestcell)+1} = ybest;
        
    % Compute aggregate sparse linear model of dominant signal based on
    % data in ybest{K+1}
        % Initialize sparse linear model
    sumk = [];
        % Build K components of sparse linear model
    for k = 1:K
        t = datacell{k}';
        sumk = [sumk t*ybest{k}];
    end
    
        % Compute sparse linear model
    sumk = sum(sumk,2)';
    
        % Accumulate sparse linear models
    sumkmat = [sumkmat;sumk];
    
    % Plot sparse linear model with dominant common signal 
    figure; 
    for k = 1:K+1

        if k < K+1

            if nrowvec(k) == 1
                subplot(K+1,1,k); plot(ybest{k},'ro'); grid on; %ylim([0 0.001]);
            else
                subplot(K+1,1,k); plot(ybest{k},'.'); grid on; %xlim([0 nrowvec(k)+1]);
            end

        else
            subplot(K+1,1,K+1); plot(sumk,'.-'); grid on; %xlim([0 N+1]);
            %hold on; plot(ybest{K+1},'or'); 
        end

    end
    
    % Assemble signifcant indices, gene list and signature matrices
    sindx = cell(K,1);
    siggenes = cell(K,1);
    siglabels = cell(K,1);
        
    for k = 1:K
        indxk = find(abs(ybest{k}));
        sindx{k} = indxk;
        siggenes{k} = datacell{k}(indxk,:);
        siglabels{k} = rownames{k}(indxk);
    end
    
    sigindxcell{length(sigindxcell)+1} = sindx; %cell2mat(sindx);
    sigmatcell{length(sigmatcell)+1} = siggenes; %cell2mat(siggenes);
    siglabelscell{length(siglabelscell)+1} = siglabels;
    
    % Do clustergrams of each realized signature matrix
    for k = 1:K
        p = sz(1);
        if p > 1
            clustergram(siggenes{k});
        end
    end   
 
    % Display FDR table after selection of sparse solution
    fdrtable
    
    % Prompt for another solution from FDR table
    prompt1 = 'Select another row # from FDR table (Enter -1 to exit): ';
    % Store optimal line from FDR table
    result = input(prompt1);
    % Exit if while-loop is prompt1 < 0
    
    if result < 0
        pflag = 0;
    end

end

% Build output structure y
y=struct('nrowvec',nrowvec,'fdrmat',fdrmat,'nsig0mat',nsig0mat,...
    'fdrtab',fdrtab,'fdrtable',fdrtable,'fdrbestmat',fdrbestmat,...
    'ybestcell',{ybestcell},'sigmatcell',{sigmatcell},...
    'sigindxcell',{sigindxcell},'siglabelscell',{siglabelscell},...
    'sumkmat',sumkmat);

disp(['Exiting JAMMIT with ',num2str(count),' models selected!!']);

return