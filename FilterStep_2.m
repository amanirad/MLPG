%
% Implements the filter step
%
% input RBF weights, RBF centers, RBF covariance matrices (eps), auxdata, and data covariance 
function [n_weights,n_m, eps_c_n ,pYY] = FilterStep_2(lambda, m, eps_c, data, Gamma)

M = length(m);
%K = length(auxdata);
%t = auxdata{1};
%data = auxdata{2};

% local mean and covariance for each RBF center
% 1D for now
n_m  = zeros(M,1);
n_C  = zeros(M,1);
n_weights = zeros(M,1);


% 
% Convert the RBF base to a weighted normal distribution.
%
%

% C = covariance
C     =  1./(  2*eps_c.^2  );
alpha_tilde = lambda./eps_c.*sqrt(pi); % converted weight see notes 

check_int_in=sum(alpha_tilde)
max(alpha_tilde)
min(alpha_tilde)
% RBF_NomralDistribution_Relation

% Disabled
% Iterate over data points
%for k=1:1 %K
k=1;

    %loop over RBFs
for i=1:M
    % weight ~ N( m_1, m_2, C_1 + C_2) 
    % Note that m_1 = data, m_2 = RBF center
    % Note that normpdf requires standard deviation input, sigma.
    % and not covariance.
    %
    % NORMPDF(X,MU,SIGMA), mu = mean, sigma =std dev
    % Assumption H=1 in y=Hx + random. 
    %
    %
    % Eq 18a where alpha in paper == alpha_tilde here. 
    %
    % Eq 18a, m_{n-1|n-2} = RBF center = m(i) here 
    % Gamma = in 18a = Gamma(k) here
    n_weights(i) = normpdf(data(k),m(i), sqrt( Gamma(k) + C(i) )  )*alpha_tilde(i);  

    % New covariance value ( C_1^-1 + C_2^-1 )^-1
    C_d_i = inv(Gamma(k));
    C_i = inv(C(i));
    n_C(i) = inv( C_d_i  + C_i );
        %
    n_m(i) = n_C(i)*( C_d_i*data(k)  + C_i*m(i)  );



end

% Calculate RBF epsilon weight from (Co) variance
eps_c_n = 1./sqrt(2*n_C);

pYY = sum(n_weights);
mma=max(n_weights)
mmi=min(n_weights)

n_weights = n_weights.*eps_c_n/sqrt(pi)./sum(n_weights);%
test=sum(n_weights./eps_c_n)*sqrt(pi)
%lambda_n(i)./sum(lambda_n).*scale; 
% Integral int_{-\infty}^\infty exp(-eps^2 (x- c)^2 ) dx = sqrt(pi) / eps
% thus we scale back with eps/sqrt(pi) 

%n_weights = eps_c'./sqrt(pi).*n_weights./sum(n_weights)

%end




