%
% Generate CIR distribuion from noncentral chi-square distribution.
%
% dX = a(b - X)dt + sigma*sqrt(X) dW
%
% The CIR pdf is related to the Chi distribution with 2q+2 degrees
% of freedom and non-centrality parameter 
% Usage: 
% Density = cirpdf(X_t , t,  X_{t-1},t-1, a, b, sigma);
% where (X_t, t) is the evaluation point at the final time
% given that the process started in (X_{t-1}, t-1)
% 
function res = cirpdf(y,t, x,s, a,b,sigma)


d = exp(-a*(t-s));
c = 2*a/(sigma.^2.*(1-d));
q = 2*a*b/(sigma.^2)-1;

% Transformed variable:
z = 2*c*y;

% Non centrality:
lambda = 2*c*x*d;

% Degrees of freedom:
dg = 2*q+2;


% Check feller condition
%feller = 2*a*b-sigma.^2;
%if(feller<0)
if(not(2*a*b>sigma.^2))    
%    warning('Feller condition not fullfilled')
res = 1e-100; % Return a very small probability
else

res = 2*c*ncx2pdf(z,dg,lambda);

end
