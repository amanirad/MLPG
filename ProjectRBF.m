%
% By Josef H??k. 
% Projects RBF 2 on RBF 1 and returns weights for RBF 1
% 
% Galerkin projection is used: 
%
%   < lambda_i varphi(x_c^i, eps_i) , lambda_j varphi(x_c^j, eps_j >
%   lambda_new  = lambda_i varphi(x_c^i, eps_i) 
%  
%
% Input: RBF 1 params,  RBF 2 params.
%    x_center, epsilon,  weight2, x_center2, epislon2 
% Output: 
% new weights on the RBF 1 space
%  
%
function [lambda] = ProjectRBF( x_c1, eps_c1, w2, x_c2, eps_c2,foo)



N_c1 = length(x_c1);
N_c2 = length(x_c2);

% Fit pdf function to a sum of RBFs using the Galerkin approach from my technical report
% and the results in Appendix A in Elisabeths and Katarinas paper.
%
% Build RBF gramian 
Gram = zeros(N_c1, N_c2);
rhs = zeros(N_c2, 1);

%
% Note eps_c1 = eps and eps_c2 = alpha 
% in the documentation
%

for i = 1:N_c1
    for j = 1:N_c2
        
%          eps_scaled = eps_c1(i)*eps_c1(j)/( eps_c1(i)^2 + eps_c1(j)^2 );
%        Gram(i,j) = (2*eps_scaled)^(1/2)*exp(-eps_scaled*eps_c1(i)*eps_c1(j).*(x_c1(i) - x_c1(j)).^2 );

        
       eps_scaled = eps_c1(i).^2*eps_c1(j).^2/( eps_c1(i)^2 + eps_c1(j)^2 );
        Gram(i,j) = (  pi./( eps_c1(i).^2 +   eps_c1(j).^2 ) ).^(1/2)*exp(-eps_scaled.*(x_c1(i) - x_c1(j)).^2 );
%        Gram(i,j) = (  2*eps_scaled).^(1/2)*exp(-eps_scaled.*(x_c1(i) - x_c1(j)).^2 );

        % Analytical calculation of the right hand side
        eps_tmp = eps_c1(i).^2*eps_c2(j).^2/( eps_c1(i)^2 + eps_c2(j)^2 );
        rhs_tmp = (  pi./( eps_c1(i).^2 +   eps_c2(j).^2 ) ).^(1/2)*exp(-eps_scaled.*(x_c1(i) - x_c2(j)).^2 );

        
        % SLASK
        %        rhs_tmp = ( 2*eps_scaled).^(1/2)*exp(-eps_scaled.*(x_c1(i) - x_c2(j)).^2 );
        %        exp(-eps_c1(i)^2.*(x - x_c1(i)).^2 );
    rhs(i) =  rhs(i) +  w2(j)*rhs_tmp;%trapz(x, foo.*rbf);

    end
    
%
% Build RHS
% the Galerking projection of pdf onto the RBFS
% Here we just use numerical trapz integration     
% It should be over index j here but we assume N_c1 = N_c2 so it does not
% matter for the moment


  %      eps_scaled = eps_c2(i).^2*eps_c2.^2./( eps_c2(i)^2 + eps_c2.^2 );
 %       rhs_tmp= w2(i).*(  pi./( eps_c2(i).^2 +   eps_c2.^2 ) ).^(1/2).*exp(-eps_scaled.*(x_c2(i) - x_c2).^2 );
%
%       rhs(i) = sum(rhs_tmp);
%rhs(i) = sqrt(2*pi)./(2*eps_c2(i)).*sum(w2);
    
    
% DISABLED TEST
% Build RHS
% the Galerking projection of pdf onto the RBFS
% Here we just use numerical trapz integration   
%x = linspace(0,4,100);
%rbf = exp(-eps_c1(i)^2.*(x - x_c1(i)).^2 );
%rhs(i) = trapz(x, foo.*rbf);
    
end
%whos    



% Calculate the new RBF weights with pseudo inverse
lambda = pinv(Gram)*rhs;

% Rescale so that the PDF integrates to unity
%  scale =;
  lambda = lambda.*eps_c1/sqrt(pi)./sum(lambda);

