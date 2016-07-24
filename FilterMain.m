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
load CIRDataSet1000
% we now have t and data arrays


FIG = 0

%[lambda, x_c, eps_c, DATA] = INIT()

%
% Time stuff
%
N_data = length(data)
t1 = diff(t);
t1 = t1(1); % assuming regular grid
t1=10;
t0 = 0;



%
% Define parameters
% Data params: a=5, b = 0.05, sigma = 0.15
a = 0.5
b = 0.05
sigma = 0.15
model_params =[a,b,sigma]

% X max for spatial grid used by DummyRBF
xmax = 2*max(data)

%
% Now define RBF node centers
%
N_c = 64;
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
pdf_state = zeros(length(xls), N_data);

%
% Initialize matrices for the solver step
%
phi = 'gs';
[Ael,Aem,Abl,Abm,Bel,Bem]=CIRMatrices(xi,xb,xls,phi,eps_c,sigma,a,b); 
M = 100; % The number of time step
[k,beta0,beta1,beta2]=BDF2coeffs(t1,M); 
S=BDF2Matrices(beta0,Ael,Aem,Abl,Abm,Bel,Bem);
%
% Initial condition, a Gaussian around the first data point
%
sigma_d = sqrt(0.0008);
u0 = normpdf(xls,x0,sigma_d);
if(FIG)
figure(5)
plot(xls,u0), hold on
end
%
% Main loop over data
%
fprintf('Working......')
for i=2:N_data
   %
   % Call RBF Fokker-Planck solver
   %
   BCrhs = @(t) CIRrhs(xb,t1,eps); % Note time is meaningless now, here 
   [lambda,mu]=BDF2LSSolver(S,k,beta1,beta2,u0,BCrhs);
   %   [lambda, x_c, eps_c] = DummyRBF( x_c, eps_c, t1, t0, x0,xmax,
   %   model_params);
   lambda1=[mu(1); lambda; mu(2)];
   rxc = xcdist(xls,x_c);
   % Kolla om eps_c_n rad eller kolumnvektor och om det funkar i RBFmat.
   A = RBFmat(phi,eps_c,rxc,'0');
   u1 = A*lambda1;
   if(FIG)
   figure(6),plot(xls,u1,'b'),hold on
   end
   %
   % read data point
   %
   test_data = data(i);
   cov_d = sigma_d^2;

   %
   % Call Filter update step
   % 
   lambda=[mu(1); lambda; mu(2)];
   [lambda_n, x_c_n, eps_c_n] = FilterStep(lambda, x_c, eps_c*ones(size(x_c)), test_data, cov_d);

   
   %
   % Project weights back to the original RBF grid with original x_c, eps_c
   %
   %lambda_foo= ProjectRBF( x_c, eps_c, lambda_n, x_c_n, eps_c_n) ;
   % Lambda foo should go back into DummyRBF as initial condition
   rxc = xcdist(xls,x_c_n);
   % Kolla om eps_c_n rad eller kolumnvektor och om det funkar i RBFmat.
   A = RBFmat(phi,eps_c_n,rxc,'0');
   u0 = A*lambda_n;
   
   
   PrintSpinBar(i, N_data) 
   
   %
   % Calculate predicted states
   %
   %
   % x = true state (unknown)
   % x_hat = E[x]  =  predicted state from the output from the Fokker-Planck
   % solver 
   %
   % x_tilde = E[x_hat] = predicted state from measurement update
   
   % approximate calulation of the expected value
   x_hat(i) = trapz(xls, xls.*u1);
   x_tilde(i) = trapz(xls, xls.*u0);   
   
   
   %
   %  Error estimates
   %
   %
   % Errors
   % E[( x_tilde - E[x_hat] )^2]
   % 
   %
   Err(i) = mean( (x_hat - x_tilde).^2 ) ;
   
   
   %
   % Plot the filtered probability
   %
   if(FIG)
   figure(5), 
    plot(xls,u0,'r');
    plot(xls,normpdf(xls,test_data,sigma_d),'--')
   end
   % Save PDF  
   pdf_state(:,i) = u0;  
end



%
% Surfplot the PDFs for the observations
%
figure
surf(pdf_state)
grid on
box on
shading interp

figure(44),
plot(data, '.')
hold on
plot(x_hat, 'ro')
plot(x_tilde, 'xc')
legend('Data', 'Predicted state before measurement update', 'Predicted state after measurement update')