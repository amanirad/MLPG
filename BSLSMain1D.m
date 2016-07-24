%function [xe,u,v,error,time,timerror]=BSLSMain1D(phi,ep,M,N,rg_c,Ne,rg_e,col)
function [maxnrm,finnrm,xs,time,error,timerror,ao,mem]=BSLSMain1D(phi,ep,M,N,rg_c,Ne,rg_e,col,mode,method,show)
%
% Main program for testing only the least squares part of the problem.  
%
% phi is the RBF 'mq', 'gs',...
% ep is the (constant) shape parameter
%  
% M is the number of time steps
% N is the number of center points
% rg_c(1:2) is the range over which the center points are uniformly 
% distributed
% Ne is the number of least squares points.
% rg_e(1:2) is the range for the least squares points  
%
% Some problem parameters that we do not change very often.  
% col is the color to use for the error plots.
%  
T = 1; % Time of maturity 
K = 1; % Strike price
ir = 0.05; % Risk free interest rate
sigma = 0.3; % Volatility
%
% Possible to choose color
%
if (nargin < 8)
  col='b';
end  
%
% The center points, evaluation points, and boundary points. 
%
if (strcmp(mode,'uni'))
  xc = linspace(rg_c(1),rg_c(2),N)';
elseif (strcmp(mode,'cheb'))
  xc = ChebPts(rg_c,N);
else
  error('Incorrect value of mode (node distribution)')
end 
%
% Put the boundary points or the extreme points last.
% 
xc = [xc(2:end-1);xc(1);xc(end)];
%
% If Ne=N-2; special case, we want collocation at nodes
%
if (strcmp(mode,'cheb'))
    xe = ChebPts(rg_c,Ne);
else %chebuni or uni
  xe = linspace(rg_e(1),rg_e(end),Ne)';
end

if (Ne==N-2)
  xe = xc(2:end-1);
end 

xb = [0;4*K];
%
% Always evaluate at a fine enough grid
%
xs = linspace(0,4*K,5*N)';
%  
% BDF2-coefficients to have constant matrices. 
%  
[k,beta0,beta1,beta2]=BDF2coeffs(T,M);
%
% Initiate all matrices or matrix factors that we need for time-stepping.
%
if (strcmp(method,'dir'))
  [Q,R,L,U,P,Cem,Abl,Abm,Ael,Aem,Esl,Esm]= ...
      BSLeastSquaresMatrices(ir,sigma,beta0,phi,ep,xc,xb,xe,xs);
elseif (strcmp(method,'qr'))
  [Q,R,L,U,P,Cem,Abl,Abm,Ael,Aem,Esl,Esm]= ...
      BSLeastSquaresMatrices_QR(ir,sigma,beta0,phi,ep,xc,xb,xe,xs);
else
  error('The method argument should be dir or qr')  
end    

%
% For the first two time-steps, lambda_1 and lambda_2 have special meaning.
%
mu_1 = 0;
lambda_2 = 0;
mu_2 = 0;
lambda_1 = BSInitial(xe,K);
%
%% Try alternative approach with LS-approx of initial data to get first lambda
%% Need another BDF2rhs function for this case, no need to go back to u
%% after
%
%f = BSInitial(xe,K);
%g = BSInitial(xb,K);
%f = f - Aem*inv(Abm)*g;
%S = Ael-Aem*inv(Abm)*Abl;
%lambda = S\f;
%mu = inv(Abm)*(g-Abl*lambda);
%lambda_1 = EvalV(Ael,Aem,lambda,mu);
%
% The time-stepping loop
%
for step = 1:M
  time(step) = sum(k(1:step));
  %
  % The right-hand sides
  %
  f = BDF2rhs(step,beta1,beta2,Ael,Aem,lambda_1,lambda_2,mu_1,mu_2);
  step
  % save dumpLS 
  %pause
  g = BSrhs(K,ir,xb,time(step));
  %
  % Move values from previous step
  %
  lambda_2 = lambda_1;
  mu_2 = mu_1;
  [lambda_1,mu_1,res]=LSSolve(Q,R,L,U,P,Cem,Abl,Abm,Ael,Aem,f,g);
  %
  % Compute the maximum error as a function of time.
  %
  u = Exact1D(xs,time(step),sigma,ir);
  v = EvalV(Esl,Esm,lambda_1,mu_1);
  timerror(step) = max(abs(u-v));
end
%
% The exact solution
%
u = Exact1D(xs,T,sigma,ir);
%
% Evaluate the RBF solution (to start with at the least squares points)
%
v = EvalV(Esl,Esm,lambda_1,mu_1);
error = u-v;
%
% Compute both maximum norm and financial norm
%
maxnrm = max(abs(error));
pos = find(xe >= 1/3 & xe <= 5/3);
finnrm = max(abs(error(pos)));

Nb = 2;
[ao,mem]=LSMLcosts(M,N,Nb*ones(size(N)),Ne);
%
% Plot the errors. The figures are not cleared, in case overlay is desired.
% 
if (~strcmp(show,'no'))
figure(1)
H=plot(xs,error,col);
%H=semilogy(xe,abs(error),col);
set(gca,'FontSize',18,'LineWidth',2)
set(H,'LineWidth',2)
xlabel('x')
ylabel('|Error|')

figure(2)
H=plot(time,timerror,col);
%H=semilogy(time,timerror,col);
set(gca,'FontSize',18,'LineWidth',2)
set(H,'LineWidth',2)
xlabel('time')
ylabel('Max error')
end