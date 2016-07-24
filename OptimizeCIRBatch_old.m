%
%
% Optimize Likelihood estimate
%



close all
clear all

parfor nsim=1:100

rand('state', 10)
randn('state',322)


OpPath = {'~/Code/Subroutines/','~/Code/RBFsupport/basisfunctions/','~/Code/RBFsupport/interpolation/', ...
	  '~/Code/RBFsupport/rbfqr/','~/Code/TestScripts/','~/Code/Optim/', '~/Code/FigScripts/','~/Code/'}

addpath(OpPath{:})
load InvCIRBigDataSet1000.mat
AllData=data;
AllT=t;



%
% Find optimal parameters
%
ic_val = [10 0.01 0.1]
% Xact values:

%%
thetaSave=zeros(100,3);
timeSave=zeros(100,1);
% a = 5
% b = 0.05
% sigma = 0.15


    data=1./AllData(nsim,101:end);
    t=AllT(101:end);
    x0  = data(1); %Start position
    % X max for spatial grid used by DummyRBF
    xmax = 2*max(data)
    %
    % Time stuff
    %
    N_data = length(data)
    t1 = diff(t);
    t1 = t1(1); % assuming regular grid
    t0 = 0;
    
    
    
    %
    % Define RBF node centers
    %
    N_c = 256;
    x_c  = linspace(0, xmax, N_c)';
    %
    % Choose how much the RBFs should decay between two nodes
    %
    s = 0.8;
    %
    % Compute the corresponding epsilon value
    %
    d = xmax./(N_c-1);
    eps_c=2*sqrt(-log(s))./d;
    
    %
    % Separate interior nodes xi, boundary nodes xb, and define least squares points
    %
    xi = x_c(2:end-1);
    xb = x_c([1 end]);
    
    xls = linspace(0,xmax,4*N_c)';
    %
    % DEfine a grid used for plotting
    %
    N = 100
    x= linspace(0,xmax,N);
    
    
    
    t00 = tic ;
    myop = optimset('MaxIter', 1, 'Diagnostics', 'on', 'Display', 'iter');
   [OptVals, ~, ~, out1] = fminsearch(  @(theta) CalculatePDF(theta, xi, xb, xls, t1, x0, x_c, data, N_data, eps_c), ic_val, myop)
    t0= toc(t00);
    
    
    thetaSave(nsim,:)=OptVals;
    timeSave(nsim)=t0;
    OutData = ['CIROptNode', num2str(nsim), '.mat'];
    
    save(OutData)
    
end

