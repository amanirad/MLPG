%
%
% Input theta xi, xb, xls,  eps_c, t1, x0,x_c,data, N_data where theta contains the parameters a,b, sigma
%
% Returns maximum of PDF
function PDF= CalculatePDF(theta, xi, xb, xls, t1, x0, x_c, data, N_data, eps_c)

a = theta(1)
b = theta(2)
sigma = theta(3)


%
% If feller condition is broken simply return a large number
%
if(2*a*b<sigma)
    PDF  = 1e10;
    return;
end
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
%sigma_d = sqrt(0.0008);
sigma_d = sqrt(0.00001);
u0 = normpdf(xls,x0,sigma_d);


PDF = 0;
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
   %
   % Output: P(X_n | Y_{1:n-1}), f?r initial condition, P(X_{n-1}| Y_{1:n-1})
   % ekvation 8, 9 
   %

   
   
   %
   % read data point
   %
   test_data = data(i);
   y_n = test_data;
   cov_d = sigma_d^2;

   %
   % Call Filter update step
   % 
   lambda=[mu(1); lambda; mu(2)];
   
   
   %
   % Plot some figures
   %
      rxc = xcdist(xls,x_c);
   A = RBFmat(phi,eps_c,rxc,'0');
   u1 = A*lambda; 
%   if(i==N_data)
%    figure(44), plot(xls, u1,'r')
%    hold off
%   end
%hold on
   
   [lambda_n, x_c_n, eps_c_n] = FilterStep(lambda, x_c, eps_c*ones(size(x_c)), test_data, cov_d);
   %
   % Output: P(X_{n-1}| Y_{1:n-1}) 
   % ekvation 10, 11 d?r jag har anv?nt ekvation 17, 18a-c
   %
   
   
   
   %
   % Project weights back to the original RBF grid with original x_c, eps_c
   %
   rxc = xcdist(xls,x_c_n);
   % Kolla om eps_c_n rad eller kolumnvektor och om det funkar i RBFmat.
   A = RBFmat(phi,eps_c_n,rxc,'0');
   u0 = A*lambda_n; 
%     if(i==N_data)
%    figure(45), plot(u0)
%    hold off
%     end
%isreal(u0)

   % Evaluate the probability density at the data point
 %  rx_data = xcdist(test_data, x_c_n);
%   A_data = RBFmat(phi, eps_c_n, rx_data, '0');
%   pdf_data = A_data*lambda_n;
   

    %
    % Log likelihood: l(theta) = sum_i log P(Y_n | Y_{1:n-1})
    %
    %
    %  P(Y_n | Y_{1:n-1}) = int P(Y_n| X_n)P(X_n | Y_{1:n-1}) dX_n
    %
    %  = sum_j \tilde{\gamma}_j int (...)
    %
    % Y_n = H*X_n + dW =>
    % P(Y_n| X_n) = normpdf(y_n, H*x_n, sigma_d) 
    % H*x_n mean and sigma_d standard deviation evaluated at y_n
    %
    % P(Y_n |Y_{1:n-1}) = \sum_i \tilde{\gamma}_i int_\mathbb{R} N(y_n,
    % H*x_n, sigma_d)*\varphi(X_{n}, m_n, C_n) dX_n
    %
    % WLOG: H = 1
    %
    % This integral,
    %
    % int_\mathbb{R} N(y_n,X_n, sigma_d)*\varphi(X_{n}, m_n, C_n) dX_n
    %
    % have an analytical solution: 
    % 
    % alpha^2 = 1/(2*sigma_d^2))
    % eps = C_n
    %
    % N(y_n,X_n, sigma_d) = sqrt(1/(2*pi*sigma_d^2))* \varphi(...)
    %
    % int varphi_eps * varphi_alpha dx = 
    %
    % sqrt(1/(2*pi*sigma_d^2))*sqrt( pi/( eps^2 + alpha^2) ) 
    % .*exp( -eps^2*alpha^2/( eps^2 + alpha^2)*( y_n - m_n)^2 ))
    %
    % P(Y_n |Y_{1:n-1}) = sum_i Gamma_i* exp( )
    %
    % Gamma _i = \tilde{\gamma}_i sqrt(1/(2*pi*sigma_d^2))* sqrt( pi/( C_n^2 +  1/(2*sigma_d^2)) )    )
    %
    
    alpha2 = 1/(2*sigma_d^2);
    

    Gamma  = lambda.*sqrt(1./(2*pi*sigma_d^2)).*sqrt( pi./( eps_c.^2 +  1/(2*sigma_d^2)) );
    exp_t = exp( -eps_c.^2*alpha2./( eps_c.^2 + alpha2).*( y_n - x_c).^2 );

    
   %    rxc = xcdist(y_n,x_c);
   % Kolla om eps_c_n rad eller kolumnvektor och om det funkar i RBFmat.
   %       A = RBFmat(phi,eps_c,rxc,'0');
   %    u11 = A*Gamma; %lambda_n.*sqrt(1./(2*pi*sigma_d^2));  % Something is wrong with the weights!

   % Create a bonafide PDF, which is a function that is both positive 
   % definit and integrates to unity.
%Lambda_min =    min(lambda)
%figure(54)
%plot(lambda)
%SUM_LAMBDA = sum(lambda)

%INTEGRAL_LAMBDA_U1 = trapz(xls, u1)

%%%u1 = bonafide(xls, u1)    
%    figure(44)
%    plot(xls, u1, 'r--')
    
    % Evalutate P(X_n | Y_{1:n-1})
%%%    u1_interps = interp1( xls, u1, y_n);
    
    %
    % Note that the first value in u1 can sometimes becomes negative 
    % To resolve negative densities we lift the entire distribution 
    % to create a bonafide density
    % 
    %
    %
%    PDF = PDF  - log(u1_interps ) 
    PDF= PDF  - log( max(Gamma'*exp_t,1e-10) );

    
%    PDF = PDF - log(max(max(u0),1e-10));
  
   PrintSpinBar(i, N_data) 
  
end
PDF = PDF/N_data
% Return - PDF peak
%PDF = -max(u0);
fprintf('Likelihood: %d\n', PDF)
end
