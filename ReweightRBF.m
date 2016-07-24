%
% Least squares reweight of the RBFs
%
function lambda = ReweightRBF( x_c1, eps_c1, w2, x_c2, eps_c2)

% Problem
% | f_ij * lambda_j  - g_i |


N_c1 = length(x_c1);
N_c2 = length(x_c2);
% Form normal equation
% Start by evaluation RBF2 at RBF 1 node points
X = zeros(N_c1, N_c2);
b = zeros(N_c2, 1);
for i = 1:N_c2
    for j = 1:N_c1        
        % f_ij
        X(i,j) = exp( -eps_c1(i).^2.*(x_c1(i) - x_c2(j)).^2 );
    end
% g_i    
b(i)  = w2(i).*exp(-eps_c2(i).^2.*(x_c1(i) - x_c2(i)).^2 );
end


lambda = (X'*X)\X'*b





