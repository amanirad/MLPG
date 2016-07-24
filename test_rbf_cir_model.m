%
% Compares a RBF solution against exact probability density for the CIR model 
% By Josef H??k and Elisabeth Larsson, 2014
%
% The Cox-Ingersoll-Ross model:
%
% dX = a(b - X)dt + sigma*sqrt(X) dW
%
% The corresponding Fokker-Planck equation
%
% d_dt f(x,t) = - d/dx ( A(x) f ) + 1/2 d^2/dx^2 ( S(x)^2 f  ) 
%
% which expands to given ( A(x) = a(b-x) and S(x) = sigma*sqrt(x) ) 
% 
% d_dt f = a f + ( a*(x-b) + sigma^2 ) d_dx f + sigma^2 x /2 d^2_dx^2 f
%
% with initial condition f(x,0) = delta(x- X)  
%
% Compare against BS model
%
% d_dt f = -r f + rx  d_dx f + sigma^2 x^2 /2 d^2_dx^2 f
%

%Params:
a = 2
b = 4
sigma = 1
%Grid
x = linspace(0, 10, 100);

% Density = cirpdf(X_t , t,  X_{t-1},t-1, a, b, sigma);
dens = cirpdf(x, 0.1, 2, 0, a, b, sigma);



figure(1) 
plot(x,dens)