%function [xs,u,v,error,time,timerror]=BSLSMLMain1D(phi,ep,M,N,rg_c,Ne,rg_e,col)
function [inrm,maxnrm,finnrm,xs,time,error,timerror,ao,mem]=...
      BSLSMLMain1D(phi,ep,M,N,rg_c,Ne,rg_e,col,mode,method,show)

%
% Main program for testing least squares together with multilevel.
% With the hope of killing of some oscillations to the left  
%
% phi is the RBF 'mq', 'gs',...
% ep(1:L) is the (constant) shape parameter for each grid
%  
% M is the number of time steps
% N(1:L) is the number of center points for each grid
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
L = length(N);
for ell=1:L
  if (strcmp(mode,'uni'))
    xc{ell} = linspace(rg_c(1),rg_c(2),N(ell))';
  elseif (strcmp(mode,'cheb'))
    xc{ell} = ChebPts(rg_c,N(ell));
  else
    error('Incorrect value of mode (node distribution)')
  end
   
  %
  % Put the boundary points or the extreme points last.
  % 
  xc{ell} = [xc{ell}(2:end-1);xc{ell}(1);xc{ell}(end)];
end
xe = linspace(rg_e(1),rg_e(end),Ne)';
if (length(N)==1 & Ne==N(1)-2)
  % To make sure we can get the collocation result
  xe = linspace(rg_e(1),rg_e(end),N(1))';
  %xe = xe(2:end-1);
  xe = xe(2:end);
end 
%
% In 1D, there are just two boundary points, simplifying matters.
%
Nb = 2;
xb = [0;4*K];
%
% Always evaluate at a fine enough grid
%
xs = linspace(0,4*K,5*N(L))';
%  
% BDF2-coefficients to have constant matrices. 
%  
[k,beta0,beta1,beta2]=BDF2coeffs(T,M);
%
% Initiate all matrices or matrix factors that we need for time-stepping.
%
for ell=1:L
  if (strcmp(method,'dir'))    
    [Q{ell},R{ell},Lf{ell},Uf{ell},P{ell},...
     Cem{ell},Abl{ell},Abm{ell},Ael{ell},Aem{ell},Esl{ell},Esm{ell}]= ...
        BSLeastSquaresMatrices(ir,sigma,beta0,phi,ep(ell),xc{ell},xb,xe, ...
                               xs);
  elseif (strcmp(method,'qr'))
    [Q{ell},R{ell},Lf{ell},Uf{ell},P{ell},...
     Cem{ell},Abl{ell},Abm{ell},Ael{ell},Aem{ell},Esl{ell},Esm{ell}]= ...
        BSLeastSquaresMatrices(ir,sigma,beta0,phi,ep(ell),xc{ell},xb,xe, ...
                               xs);
  else
    error('The method argument should be dir or qr')      
  end
end
xe
%
% For the first two time-steps, lambda_1 and lambda_2 have special meaning.
%
for ell=1:L
  mu_1{ell} = 0;
  lambda_2{ell} = 0;
  mu_2{ell} = 0;
  lambda_1{ell} = zeros(Ne,1);
end  
lambda_1{L} = BSInitial(xe,K);
f = zeros(Ne,1);
%
% The time-stepping loop
%
for step = 1:M
  figure(4),clf
  time(step) = sum(k(1:step));
  %
  % The right-hand sides
  %
  f(:) = 0;
  for ell=1:L
    f = f + BDF2rhs(step,beta1,beta2,Ael{ell},Aem{ell}, ...
                    lambda_1{ell},lambda_2{ell},mu_1{ell},mu_2{ell});
  end  
  g = BSrhs(K,ir,xb,time(step));
  %
  % Move values from previous step
  %
  lambda_2 = lambda_1;
  mu_2 = mu_1;
  %
  % First solve on the coarsest grid. f and g are the global ones.
  %
  [lambda_1{1},mu_1{1},res]=LSSolve(Q{1},R{1},Lf{1},Uf{1},P{1}, ...
                              Cem{1},Abl{1},Abm{1},Ael{1},Aem{1},f,g);
  %
  % For the rest of the grids, f is the residual and g is zero (in 1D).
  %
  semilogy(abs(res)),hold on
  res_l(1,step)=max(abs(res));
  for ell=2:L
    f = res;
    g(:) = 0;
    [lambda_1{ell},mu_1{ell},res]= ...
        LSSolve(Q{ell},R{ell},Lf{ell},Uf{ell},P{ell}, ...
                Cem{ell},Abl{ell},Abm{ell},Ael{ell},Aem{ell},f,g);
    res_l(ell,step)=max(abs(res));
    semilogy(abs(res)),hold on
  end  
  %
  % Compute the maximum error as a function of time.
  %
  u = Exact1D(xs,time(step),sigma,ir);
  v = EvalV(Esl{1},Esm{1},lambda_1{1},mu_1{1});
  %figure(3),clf
  %plot(xe,v),hold on
  for ell=2:L
    err_l(ell-1,step) = max(abs(u-v));
    semilogy(abs(u-v),'r--'),hold on
    v0=v;
    v = v + EvalV(Esl{ell},Esm{ell},lambda_1{ell},mu_1{ell});
    corr_l(ell-1,step) = max(abs(v-v0));
    %plot(xe,v)
  end
  semilogy(abs(u-v),'r--'),hold on
  err_l(L,step) = max(abs(u-v));
  timerror(step) = max(abs(u-v));
  %pause(0.1)
end
%
% The exact solution
%
u = Exact1D(xs,T,sigma,ir);
%
% Evaluate the RBF solution (to start with at the least squares points)
%
% Already done?
%v = EvalV(Ael,Aem,lambda_1,mu_1);

error = u-v;
%
% Compute both maximum norm and financial norm
%
maxnrm = max(abs(error));
pos = find(xs >= 1/3 & xs <= 5/3);
finnrm = max(abs(error(pos)));

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
%
% Compute the integrated timeerror
%
h = time(2:end)-time(1:end-1);
inrm=sum(h/2.*(log10(timerror(1:end-1))+log10(timerror(2:end))));

if (L>1)
  figure(3),semilogy(time,err_l,'-',time,corr_l,'--',time,res_l,'-.')
end

