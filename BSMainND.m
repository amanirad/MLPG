% Beräkna residual och fel i 1-D och 2-D och beräkna sen residualen i
% fler-d. Gör en graf. Eventuellt analytisk uppskattning av fel från
% residual. Om jag enkelt kan göra det från heat eq, så kan jag göra det.
% Skala ep med 1/D. 2-D till 3-D så gick residualen ner. Kanske till och
% med kan göra Cauchy-problemet eftersom vi inte har egentliga ränder
% efter log-transform.
%
% Se till att spara residualen. Beräkna två olika lösningar och
% jämför. Second oder?xs
function [xs,sol,time,res_l,R]=BSMainND(phi,ep,M,N,d,ctype,etype,method,show,col,Nls,Ns);
%
% Some problem parameters that we do not change very often.  
%
T = 1; % Time of maturity 
K = 1; % Strike price
ir = 0.05; % Risk free interest rate
sigma = 0.05*ones(d);
sd = diag(diag(sigma));
sigma = sigma - sd + 0.3*eye(d); % Volatility
%
% Compute nodes, collocation/least squares points, and error evaluation points 
%
L = length(N);

for ell=1:L
  %%%%%%%%%%%%%%%%%%% redefine N?
  [xi,b1,b2] = GenNDPoints(N(ell),d,ctype);
  %
  % Put the boundary points or the extreme points last.
  % 
  xc{ell} = [xi;b1;b2];
  xb{ell} = [b1;b2];
  disp(['Interior points grid ' num2str(ell) ': ' num2str(size(xi,1))])
  disp(['Boundary points grid ' num2str(ell) ': ' ...
        num2str(size(b1,1)+size(b2,1))])
end
xs = GenNDPoints(Ns,d,etype);
%
% We will always use the least squares variant. The collocation cases are
% not implemented other than through making Nls=N;
%

if (L==1 & Nls==N)
  xls = xi; %%%%%%%%%%%%%%%%% Does not work yet. Wrong number?
else
  xls = GenNDPoints(Nls,d,etype);
  Nls = size(xls,1);
end
%
% We do not have a reference solution for the d-dimensional problem. We
% could possibly make one? Ask Lina about FD-code. Or serious Monte Carlo simulation?
%
%
% BDF2-coefficients to have constant matrices. 
%  
[k,beta0,beta1,beta2]=BDF2coeffs(T,M);
%
% Initiate all matrices or matrix factors that we need for time-stepping.
%
[Q,R,Lf,Uf,P,Cem,Abl,Abm,Ael,Aem,Esl,Esm,Albml,Albmm]= ...
    BSLeastSquaresMatrices(method,ir,sigma,beta0,phi,ep,xc,xb,xls,xs);
%
% For the first two time-steps, lambda_1 and lambda_2 have special meaning.
%
for ell=1:L
  mu_1{ell} = 0;
  lambda_2{ell} = 0;
  mu_2{ell} = 0;
  lambda_1{ell} = zeros(Nls,1);
end
lambda_1{L} = BSInitial(xls,K);
f = zeros(Nls,1);%
res = 0; timeres=0;
% The time-stepping loop
%
for step = 1:M
  time(step) = sum(k(1:step));
  %
  % The right-hand sides
  %
  f(:) = 0;
  for ell=1:L
    f = f + BDF2rhs(step,beta1,beta2,Ael{ell},Aem{ell}, ...
                    lambda_1{ell},lambda_2{ell},mu_1{ell},mu_2{ell});
  end  
  g = BSrhs(K,ir,xb{1},time(step));
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
  res_l(1,step) = max(abs(res)); 
    %
  % For the rest of the grids, f is the residual and g is the difference.
  %
  for ell=2:L
    f = res;
    g = BSrhs(K,ir,xb{ell},time(step));
    for m=1:ell-1
      g = g - Albml{ell}{m}*lambda_1{m}-Albmm{ell}{m}*mu_1{m};  
    end  
    [lambda_1{ell},mu_1{ell},res]= ...
        LSSolve(Q{ell},R{ell},Lf{ell},Uf{ell},P{ell}, ...
                Cem{ell},Abl{ell},Abm{ell},Ael{ell},Aem{ell},f,g);
    control =max(abs(lambda_1{ell}))
    res_l(ell,step)=max(abs(res));
  end  
  %
end
%
% Evaluate the RBF solution at the points xs
%
sol = EvalV(Esl{1},Esm{1},lambda_1{1},mu_1{1});
for ell=2:L
  sol = sol + EvalV(Esl{ell},Esm{ell},lambda_1{ell},mu_1{ell});
end  
%
% Compute only the error
%
%error = (rsol(:)-vv); %./(max(rsol(:),1));
%rg_err = [min(abs(error(pos2))) max(abs(error(pos2))) min(abs(error(pos1)))]; 
%error = reshape(error,size(ss1));
%vv = reshape(vv,size(ss1));

%if (~strcmp(show,'no'))
  %
  % Plot the errors. The figures are not cleared.
  % 
  %  figure(1),clf
  %z = log10(abs(error));
  %rg_err = log10(rg_err);
  %cl = rg_err(1):(rg_err(2)-rg_err(1))/20:rg_err(2);
  %  cl =[rg_err(3) cl]; 
  %[C,H]=contourf(ss1,ss2,z,cl);
  %set(H,'EdgeColor','none')
  %set(gca,'FontSize',18,'LineWidth',2)
  %xlabel('x_1')
  %ylabel('x_2')
  %title('Error')
  %H=colorbar;
  %set(H,'FontSize',18,'LineWidth',2)
  %hold on
  %sv1=0:0.01:2;
  %sv2=2-sv1;
  %H=plot(sv1,sv2,'k--');
  %set(H,'LineWidth',2)
  %end