%
% Test filter step
%
% This script test and validates the FilterStep function
% By Josef H??k, josef.hook@it.uu.se
% 2014-10-30
%
%
close all
clear all


%
% load some data
%
load Sparse_CIRDataSet1000
% we now have t and data arrays


%%
%
% Part I
% Build an RBF approximate represenation of the probability density
% function from the analytical PDF for the CIR process
% using the Galerkin projection technique.
%
% This part simulates an output from the RBF Fokker-Planck solver
%
%
%

%
% Define grid and parameters
%
N = 200;
x = linspace(0, 10, N);
t = 0.1
t0 = 0
x0  = 2
a = 5
b = 0.4
sigma = 1

%
% Remember Feller condition 2*a*b > sigma
%
pdf =  cirpdf(x , t,  x0,t0, a, b, sigma);
figure(22)
plot(x, pdf)

%
% Now define RBF node centers
%
N_c = 40
x_c  = linspace(0, 10, N_c)
eps_c = 4.2*ones(N_c,1);

% eps_c is sensitive to the distance between the observed data point
% and the initial point
%
% Examples eps_c =1.2 and data =2.2 gives good result while
% data =4 breaks down
% If one increases eps_c =4.2 then data =3  give good result x0=2
%
% eps_c = 6.2 good for data =4, x0 =2
%
% NOTE:
% There seems to be an empirical relation between the relative distance
% between x0, data and the eps_c. This probably has to do with the level 
% of intersection of the Gaussian bells. Little intersection, little
% information sharing etc...
% 
%
%

%
% Fit pdf function to a sum of RBFs using the Galerkin approach from my technical report
% and the results in Appendix A in Elisabeths and Katarinas paper.
%
% Build RBF gramian
%
Gram = zeros(N_c);
rhs = zeros(N_c, 1);
for i = 1:N_c
    for j = 1:N_c
        eps_scaled = eps_c(i)*eps_c(j)/( eps_c(i)^2 + eps_c(j)^2 );
        Gram(i,j) = (pi./( eps_c(i).^2 +   eps_c(j).^2 ))^(1/2)*exp(-eps_scaled*eps_c(i)*eps_c(j).*(x_c(i) - x_c(j)).^2 );
    end
    
    %
    % Build RHS
    % the Galerking projection of pdf onto the RBFS
    % Here we just use numerical trapz integration
    %
    rbf = exp(-eps_c(i)^2.*(x - x_c(i)).^2 );
    rhs(i) = trapz(x, pdf.*rbf);
    
    
    % plot each RBF
    if(i==4)
        figure(22)
        hold on
        plot(x, rbf, 'r--')
        drawnow
    end
    
    
end

%
% Calculate the RBF weights
% pdf_approx = sum_n^N lambda_n * phi_eps( x, x_c )
%
lambda = Gram\rhs

%
% Normalize lambda such that the integral of PDF is == 1
%

%
% Build and Plot the approximate pdf
%
pdf_approx = zeros(1,N);
for i = 1:N_c
    pdf_approx = pdf_approx + lambda(i).*exp(-eps_c(i)^2.*(x - x_c(i)).^2 );
    
end
plot(x, pdf_approx, 'k*-')



%%
%
% Part II
% test the filter step
%
% We now have from part I a represtentation of the PDF in terms of a sum of
% RBFs
%
% Input is the weights lambda, the node centers x_c and the eps_c weights
% together with a data point
%
%


%
% Define some data with measure uncertainty
%
cov_d = 0.1
test_data = 3 % data(1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILTER STEP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:     auxiliary data
% Ouput:
% [n_weights,n_m, n_C ] = FilterStep(lambda, m, C, auxdata, C_data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILTER STEP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lambda_n, x_c_n, eps_c_n] = FilterStep(lambda, x_c, eps_c, test_data, cov_d)


%
%Plot data point measurement distribution
%
plot(test_data, 1, 'rx', 'linewidth', 6)
plot(x, normpdf(x, test_data, sqrt(cov_d)), 'm')

%
% Build and Plot the approximate pdf
%
pdf_filtered = zeros(1,N);
for i = 1:N_c
    
    pdf_filtered = pdf_filtered + lambda_n(i).*exp(-eps_c_n(i)^2.*(x - x_c_n(i)).^2 );
    
    if(i==4)
        rbf_n = exp(-eps_c_n(i)^2.*(x - x_c_n(i)).^2 );
        plot(x, rbf_n, '--')
    end
    
end

%
% Plot the filtered distribution
%
plot(x, pdf_filtered, 'co-')
Integral_PDF_filtered = trapz(x, pdf_filtered)


