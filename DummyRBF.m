%
% Dummy RBF solver using analytical CIR solution
% By Josef H??k, josef.hook@it.uu.se
% 
function [lambda, x_c, eps_c] = DummyRBF(  x_c, eps_c, t1, t0, x0,xmax, model_params)


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
x = linspace(0, xmax, N);
a = model_params(1);
b = model_params(2);
sigma = model_params(3);

%
% Remember Feller condition 2*a*b > sigma
%
pdf =  cirpdf(x , t1,  x0,t0, a, b, sigma);

FIG = 0;
if(FIG)
figure(22)
plot(x, pdf)
end

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
N_c=length(eps_c);
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
    
    if(FIG)

    % plot each RBF
    if(i==4)
        figure(22)
        hold on
        plot(x, rbf, 'r--')
        drawnow
    end
    end
    
    
end

%
% Calculate the RBF weights
% pdf_approx = sum_n^N lambda_n * phi_eps( x, x_c )
%
%lambda = Gram\rhs;
lambda = pinv(Gram)*rhs;
%
% Normalize lambda such that the integral of PDF is == 1
%
if(FIG)

%
% Build and Plot the approximate pdf
%
pdf_approx = zeros(1,N);
for i = 1:N_c
    pdf_approx = pdf_approx + lambda(i).*exp(-eps_c(i)^2.*(x - x_c(i)).^2 );
    
end
plot(x, pdf_approx, 'k*-')


end