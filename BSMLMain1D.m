Bfunction [inrm,maxnrm,finnrm,xe,time,error,timerror,ao,mem]= ...
      BSMLMain1D(phi,ep,M,N,rg_c,col,mode,method,show)
%
% Main program for testing multilevel alone.
%
% phi is the RBF 'mq', 'gs',...
% ep(1:L) is the (constant) shape parameter for each grid
%  
% M is the number of time steps
% N(1:L) is the number of center points for each grid
% rg_c(1:2) is the range over which the center points are uniformly 
% distributed
% col is the color to use for the error plots.
% mode is 'uni' or 'cheb' 
% 
% Some problem parameters that we do not change very often.  

% mode is uni or cheb  
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
%
% xe is here only used for evaluation purposes.
%
xe = linspace(rg_c(1),rg_c(end),5*N(L))';
%
% In 1D, there are just two boundary points, simplifying matters.
%
Nb = 2;
xb = [0;4*K];
%  
% BDF2-coefficients to have constant matrices. 
%  
[k,beta0,beta1,beta2]=BDF2coeffs(T,M);
%
% Initiate all matrices or matrix factors that we need for time-stepping.
%

if (strcmp(method,'dir'))
    [L2,U2,P2,L1,U1,P1,Cel,Cem,Abl,Abm,Ael,Aem,Esl,Esm]= ...
        BSMLMatrices(ir,sigma,beta0,phi,ep,xc,xb,xe);
elseif (strcmp(method,'qr'))
    [L2,U2,P2,L1,U1,P1,Cel,Cem,Abl,Abm,Ael,Aem,Esl,Esm]= ...
        BSMLMatrices_QR(ir,sigma,beta0,phi,ep,xc,xb,xe);
else
  error('The method argument should be dir or qr')  
end

    %
% For the first two time-steps, lambda_1 and lambda_2 have special meaning.
%
%for ell=1:L
%  mu_1{ell} = 0;
%  lambda_2{ell} = 0;
%  mu_2{ell} = 0;
%  lambda_1{ell} = zeros(Ne,1);
%end  
%lambda_1{L} = BSInitial(xe,K);
%f = zeros(Ne,1);
%
% The time-stepping loop
%
for step = 1:M
  time(step) = sum(k(1:step));
  %
  % Solve for one level at a time using collocation and explicit
  % computation of the residual
  %
  for ell = 1:L
    %
    % The right-hand sides
    %
    if (step==1)
      for m=1:L
        mu_1{m} = 0;
        lambda_1{m} = 0;
        mu_2{m} = 0;
        lambda_2{m} = 0;
      end
      lambda_1{ell} = BSInitial(xc{ell}(1:end-Nb,:),K);
    elseif (step==2)
      for m=1:L
        lambda_2{m} = 0;
        mu_2{m} = 0;
      end
      lambda_2{ell} = BSInitial(xc{ell}(1:end-Nb,:),K); 
    end 
    f = zeros(size(xc{ell},1)-Nb,1);
    for m=1:L
      f = f + BDF2rhs(step,beta1,beta2,Ael{ell,m},Aem{ell,m}, ...
                      lambda_1{m},lambda_2{m},mu_1{m},mu_2{m});
    end 
    for m=1:ell-1
      f = f - Cel{ell,m}*lambda{m} - Cem{ell,m}*mu{m};
    end  
    g = BSrhs(K,ir,xb,time(step));
    %
    % Test to see how exact conditions matter for small errors
    %
    %g = Exact1D(xb,time(step),sigma,ir);
    if (ell > 1)
      g = 0*g;
    end  
    [lambda{ell},mu{ell}]=MLSolve(L2{ell},U2{ell},P2{ell}, ...
                                  L1{ell},U1{ell},P1{ell}, ...
                                  Cem{ell,ell},Abl{ell},Abm{ell}, ...
                                  Ael{ell,ell},Aem{ell,ell},f,g);
  end  
  %
  % Move values from previous step
  %
  lambda_2 = lambda_1;
  mu_2 = mu_1;
  lambda_1 = lambda;
  mu_1 = mu;
  %
  % Compute the maximum error as a function of time.
  %
  u = Exact1D(xe,time(step),sigma,ir);
  v = EvalV(Esl{1},Esm{1},lambda_1{1},mu_1{1});
  %figure(3),clf
  %plot(xe,v),hold on
  for ell=2:L
    err1 = max(abs(u-v))
    v0=v;
    v = v + EvalV(Esl{ell},Esm{ell},lambda_1{ell},mu_1{ell});
    corr1 = max(abs(v-v0))
    %plot(xe,v)
  end  
  timerror(step) = max(abs(u-v));
  %pause(0.1)
end
%
% The exact solution
%
u = Exact1D(xe,T,sigma,ir);
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
pos = find(xe >= 1/3 & xe <= 5/3);
finnrm = max(abs(error(pos)));

[ao,mem]=CLMLcosts(M,N,Nb*ones(size(N))); % OK also for collocation


%
% Plot the errors. The figures are not cleared, in case overlay is desired.
% 
if (~strcmp(show,'no'))
figure(1)
H=plot(xe,error,col);
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

