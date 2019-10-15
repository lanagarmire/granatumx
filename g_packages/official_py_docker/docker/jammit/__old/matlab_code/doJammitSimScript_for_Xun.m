% JAMMIT simulation for a single data matrix 
%
nrows=1000;     % Number of rows of simulated data matrix
ncols=50;      % Number of columns of simulated data matrix
nsig=50;        % Signal size
nalpha=25;      % Length of alpha grid
nperms=30;      % Number of permutations for estimating FDR

% Generate noise data matrix
noisefac=1.0;
noisemat=noisefac.*randn(nrows,ncols);
% Generate signal data matrix
sigmat=zeros(nrows,ncols);
sindx=randsample([1:nrows],nsig);
% Generate simulated step signal
step=ones(1,ncols); 
step(1:ncols/2)=-ones(1,ncols/2);
sigfac=0.75;
step=sigfac.*step;
% Embedd step signal in signal data matrix
sigmat(sindx,:)=repmat(step,nsig,1);   
% Generate simulated data matrix with embedded step signal
datamat=sigmat+noisemat;
% Compute matrix SNR
snr = 10*log10(var(sigmat)/var(noisemat));
snr
% Define JAMMIT input variables
datacell={datamat};
nsuppvec=[nsig];
alpha0=0;
alpha1=max(abs(datamat(:)))+1;
alphavec=linspace(alpha0,alpha1,nalpha);
rowids={[1:nrows]};
colids={[1:ncols]};
% Run JAMMIT with FDR table
efdrout = efdr_ssvdanyN(datacell,nsuppvec,alphavec,nperms,rowids,colids);
% y=efdr_ssvdsymN({z},[nsig],[100],linspace(0,x,5),5,{[1:2000]},...
% {[1:100]},0);

% Generate plots of 1st JAMMIT solution
if ~isempty(y)
    temp=efdrout.ybestcell{:};
    u=temp{1};
    v=temp{2};
    indx=find(abs(u)>0);
    lenu=length(u);
    % Plot JAMMIT detector output
    figure; plot(u,'.'); grid on; hold on; xlim([0 lenu]);
    plot(sindx,u(sindx),'or');      % Highlight true support
    plot(indx,u(indx),'+g');        % Highlight JAMMIT support
    return;
else
    disp('Empty Solution!!');
    return;
end