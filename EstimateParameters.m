%
%
% Maximum likelihood optimization of the filter 
% 

%%%%%%%%%%%
%
% N?got ?r fel i koden den ger helt fel optimimum. 
% M?ste kolla s? att likelihooden ?r korrekt formulerad 
% Ibland f?r tex [10 0.1 0.1] [a,b, sigma] s? blir log(u) komplex vilket
% indikerar att n?got ?r galet. 
% Vilken os?kerhet ska man ha p? m?tfelet, sigma_d????
% 
%

function EstimateParameters
%
% By Josef H??k 
% 2014-10-30
% Filter Main
%
%
close all
clear all

%
% load some data
%
load Sparse_CIRDataSet1000
% we now have t and data arrays


%
% Time stuff
%
N_data = length(data)
t1 = diff(t);
t1 = t1(1); % assuming regular grid
%t1=10;
t0 = 0;



%
% Define parameters
% Data params: a=5, b = 0.05, sigma = 0.15
%a = 0.5
%b = 0.05
%sigma = 0.15
%model_params =[a,b,sigma]

% X max for spatial grid used by DummyRBF
xmax = 2*max(data)

%
% Now define RBF node centers
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

x0  = data(1) %Start position
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


%
% Find optimal parameters
%
ic_val = [10 0.01 0.1]
% Xact values:
% a = 5
% b = 0.05
% sigma = 0.15

%myop = optimoptions('fminunc','FinDiffType', 'central', 'MaxIter', 50, 'Algorithm', 'quasi-newton', 'Diagnostics', 'on', 'Display', 'iter');   
%OptVals = fminunc( @(theta) CalculatePDF(theta, xi, xb, xls, t1, x0, x_c, data, N_data, eps_c), ic_val, myop)

%'MaxIter', 50, 
myop = optimset('Diagnostics', 'on', 'Display', 'iter');
OptVals = fminsearch(  @(theta) CalculatePDF(theta, xi, xb, xls, t1, x0, x_c, data, N_data, eps_c), ic_val, myop)


end

