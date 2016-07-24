function [maxnrm,xs,error,xls,res,time,timerror,timeres,ao,mem]= ...
      CIRcompare(a_p,b_p,sigma_p,T,x0,start,ep,M,N,Ns,rg,nodetype,method,Nls)
%
% Test case: 
% CIRMain('gs',1,100,40,200,[0 20],'uni','uni',{'dir'},1,'b',80)
% Main program for all the different variants of the method 
%
% phi is the RBF 'mq', 'gs',...
% ep(1:L) is the (constant) shape parameter for each grid
%         Note that for a one grid method ep is just given as a scalar 
% M is the number of time steps
% N(1:L) is the number of center points for each grid
% Ns is the number of uniform points where the solution is evaluated
% rg(1:2) is the range over which we are computing, typiclly [0 4]
% ctype is 'cheb' or 'uni' for the distribution of node points
% etype has the same function for the points where equations are enforced
% method{1:L} is 'qr' or 'dir' for RBF-QR or RBF-Direct for each grid
%             note only the lsml version allows different
% show is 'yes' or 'no' and tells if plots of the errors should be made
% col is the color to use for the error plots.
% Nls is _optional_. If present it indicates that least squares should be used.
%-----------------------------------------
%
% Some problem parameters that we do not change very often.  
%
phi = 'gs'
% Initial time to get smooth initial condition, used by CIRinitial
t1 = start;%1e-1; 
t0 = 0; 
%x0 = 1; % initial positition for the process
ctype=nodetype;
etype=ctype;
show='yes';
col='b';

% 
% Determine if least squares are to be used
%
if (nargin<=7)
  ls = 0;
else
  ls=1;
end
%
% Compute nodes, collocation/least squares points, and error evaluation points 
% But this is only for 1-D, right?
% 
L = length(N);
for ell=1:L
  if (strcmp(ctype,'uni'))
    xc{ell} = linspace(rg(1),rg(2),N(ell))';
  elseif (strcmp(ctype,'cheb'))
    xc{ell} = ChebPts(rg,N(ell));
  elseif (strcmp(ctype,'adap'))  
    xc{ell} = AdaptivePts(rg,N(ell));
  else
    error('Incorrect value of mode (node distribution)')
  end
  %
  % Put the boundary points or the extreme points last.
  % 
  xc{ell} = [xc{ell}(2:end-1);xc{ell}(1);xc{ell}(end)];
end

if (ls)
  if (strcmp(etype,'cheb'))
    xls = ChebPts(rg,Nls);
  elseif (strcmp(etype,'uni'))
    xls = linspace(rg(1),rg(end),Nls)';
  end
  %
  % Allow a special case for Nls=N-2 with collocation and QR-factorization
  % 
  if (L==1 & Nls==N(1)-2)
    xls = xc{1}(2:end-1);
  end
else
  xls = 0;
end
%
% Boundary points
%
xb= rg(:); % For the 1-D case
Nb = 2;
%
% The evaluation points are always uniform Ns should be larger than N
%
xs = linspace(rg(1),rg(2),Ns)';
%  
% BDF2-coefficients to have constant matrices. 
%  
[k,beta0,beta1,beta2]=BDF2coeffs(T-t1,M);
%
% Initiate all matrices or matrix factors that we need for time-stepping.
%
if (ls)
  
  [Q,R,Lf,Uf,P,Cem,Abl,Abm,Ael,Aem,Esl,Esm,Albml,Albmm]= ...
        CIRLeastSquaresMatrices(method,a_p,b_p,sigma_p,beta0,phi,ep,xc,xb,xls,xs);
else
  if (strcmp(method{1},'dir'))
    [L2,U2,P2,L1,U1,P1,Cel,Cem,Abl,Abm,Ael,Aem,Esl,Esm]= ...
        CIRMLMatrices(a_p,b_p,sigma_p,beta0,phi,ep,xc,xb,xs);
  elseif (strcmp(method{1},'qr'))
      error('QR not implemented, yet') 
%    [L2,U2,P2,L1,U1,P1,Cel,Cem,Abl,Abm,Ael,Aem,Esl,Esm]= ...
 %       BSMLMatrices_QR(ir,sigma,beta0,phi,ep,xc,xb,xs);
  else
    error('The method argument should be dir or qr')  
  end     
end
%
% The time-stepping loop
%
res = 0; timeres=0;
for step = 1:M
  time(step) = t1+sum(k(1:step));
  %
  % Start with the right hand sides that do not change with the level
  %
  g = CIRrhs(xb,time(step), ep);
  %
  % Loop over the levels, solving for one at a time
  % 
  for ell = 1:L
    %
    % For the first two time-steps, lambda_1 and lambda_2 have special meaning.
    %
    if (step==1)
      for m=1:L
        mu_1{m} = 0;
        lambda_2{m} = 0;
        mu_2{m} = 0;
        lambda_1{m} = 0;
      end
      if (ls)
        lambda_1{L} = CIRInitial(xls,t1, x0, t0, a_p,b_p,sigma_p);
      else
        lambda_1{ell} = CIRInitial(xc{ell}(1:end-Nb,:),...
            t1, x0, t0, a_p,b_p,sigma_p);
      end  

    elseif (step==2 & ~ls)
      for m=1:L
        lambda_2{m} = 0;
        mu_2{m} = 0;
      end
      if (~ls)
        lambda_2{ell} = CIRInitial(xc{ell}(1:end-Nb,:),...
                t1, x0, t0, a_p,b_p,sigma_p); 
      end
    end
    %
    % The least squares rhs
    %
    if (ls & ell==1) % Is the residual in later steps
      f = zeros(Nls,1);
      for m=1:L
        f = f + BDF2rhs(step,beta1,beta2,Ael{m},Aem{m}, ...
                        lambda_1{m},lambda_2{m},mu_1{m},mu_2{m});
      end
    end
    %
    % The right-hand side for the collocation case
    % 
    if (~ls)
      f = zeros(size(xc{ell},1)-Nb,1);
      for m=1:L
        f = f + BDF2rhs(step,beta1,beta2,Ael{ell,m},Aem{ell,m}, ...
                        lambda_1{m},lambda_2{m},mu_1{m},mu_2{m});
      end 
      for m=1:ell-1
        f = f - Cel{ell,m}*lambda{m} - Cem{ell,m}*mu{m};
      end  
    end
    if (ls)
      [lambda{ell},mu{ell},res(1:Nls,ell)]= ...
          LSSolve(Q{ell},R{ell},Lf{ell},Uf{ell},P{ell}, ...
                  Cem{ell},Abl{ell},Abm{ell},Ael{ell},Aem{ell},f,g);
      timeres(step,ell) = rg(2)/(N(ell))*sum(abs(res(:,ell)));%max(abs(res(:,ell)));%sqrt(rg(2)/(N(ell)-1)*sum(res(:,ell).^2));
    else
      [lambda{ell},mu{ell}]=MLSolve(L2{ell},U2{ell},P2{ell}, ...
                                    L1{ell},U1{ell},P1{ell}, ...
                                    Cem{ell,ell},Abl{ell},Abm{ell}, ...
                                    Ael{ell,ell},Aem{ell,ell},f,g);
    end
    %
    % At all levels except the first g=0 in 1-D
    %
    g(:) = 0;
    if (ls)
      % The new right hand side is the residual
      f = res(:,ell);
    end
  end
  coeff_l=1- normcdf(0, xc{1}(1:N-Nb),  1/sqrt(2)/ep(1));
  coeff_m= 1- normcdf(0, xc{1}(N-Nb+1:N),  1/sqrt(2)/ep(1));
  pos_l=find(lambda{1}<0); pos_m=find(mu{1}<0);
  negmass(step)=sqrt(pi)/ep(1)*(sum(coeff_l(pos_l).*lambda{1}(pos_l))+...
                                sum(coeff_m(pos_m).*mu{1}(pos_m)));
  mass(step)=sqrt(pi)/ep(1)*(sum(coeff_l.*lambda{1})+sum(coeff_m.*mu{1}));
  %    [lambda{1};mu{1}]
  %    pause
  %  end
  %
  % Try Crazy approach, normalize after each step
  %
  %lambda{1}=(1/mass(step))*lambda{1};
  %mu{1}=(1/mass(step))*mu{1};
  %
  % Move the values from the previous step
  %
  lambda_2 = lambda_1;
  mu_2 = mu_1;
  lambda_1 = lambda;
  mu_1 = mu;
  %
  % Compute the maximum error as a function of time.
  %
  u = cirpdf(xs, time(step), x0, t0, a_p, b_p, sigma_p);
%  u = Exact1D(xs,time(step),sigma,ir);
  v = EvalV(Esl{1},Esm{1},lambda_1{1},mu_1{1});
  
  %figure(44)
  %plot(xs, v, 'k', xs, u, 'r--')
  %legend('RBF sol', 'Exact')
  

  
  for ell=2:L
    v = v + EvalV(Esl{ell},Esm{ell},lambda_1{ell},mu_1{ell});
  end
  timerror(step) = max(abs(u-v));
  %  if (min(lambda{1})<0 | min(mu{1})<0)
    %    figure(38),clf
    %    plot([lambda{1};mu{1}],'o')

    
end
%
% Compute the errors at the final time
%
error = u-v;
%
% Compute both maximum norm and financial norm
%
maxnrm = max(abs(error));
pos = find(xs >= 1/3 & xs <= 5/3);
finnrm = max(abs(error(pos)));
%
% Compute the integrated timeerror
%
h = time(2:end)-time(1:end-1);
inrm=sum(h/2.*(log10(timerror(1:end-1))+log10(timerror(2:end))));
%
% Compute the computationl costs
%
if (ls)
  [ao,mem]=LSMLcosts(M,N,Nb*ones(size(N)),Nls); % OK also for LS
else
  [ao,mem]=CLMLcosts(M,N,Nb*ones(size(N))); % OK also for collocation
end
%
% Plot the errors. The figures are not cleared, in case overlay is desired.
% 
if (~strcmp(show,'no'))
  figure(1)
  H=plot(xs,error,col);
  hold on
  if (ls)
    H=[H; plot(xls,res,'r--')];
  end
  %H=semilogy(xe,abs(error),col);
  set(gca,'FontSize',18,'LineWidth',2)
  set(H,'LineWidth',2)
  xlabel('x')
  ylabel('|Error|')

  figure(2)
  H=plot(time,timerror,col);
  hold on
  H=[H; plot(time,timeres,'r--')];
  %H=semilogy(time,timerror,col);
  set(gca,'FontSize',18,'LineWidth',2)
  set(H,'LineWidth',2)
  xlabel('time')
  ylabel('Max error')
  
  figure(3)
  plot(time,negmass)
  
  figure(4)
  plot(time,mass)
end

%
% Now compute the solution in the new way and compare the result.
%
[Ael,Aem,Abl,Abm,Bel,Bem]=CIRMatrices(xc{1}(1:end-2),xb,xls,phi,ep(1),sigma_p,a_p,b_p);
S=BDF2Matrices(beta0,Ael,Aem,Abl,Abm,Bel,Bem);
u0e = CIRInitial(xls,t1, x0, t0, a_p,b_p,sigma_p);
BCrhs = @(t) CIRrhs(xb,t+t1,ep);
[lambda_new,mu_new]=BDF2LSSolver(S,k,beta1,beta2,u0e,BCrhs);
diff_l = lambda{1}-lambda_new
diff_m = mu{1}-mu_new