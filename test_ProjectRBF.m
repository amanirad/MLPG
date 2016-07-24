%
% By Josef H??k
% Test the projection step as described in 
% Doc/GalerkinProjection.pdf 
%
%

% assume you have run test_FilterStep.m 


new_lambdas= ProjectRBF( x_c, eps_c, lambda_n, x_c_n, eps_c_n,pdf_filtered ) 

%
% Plot the filtered distribution in the moved RBFs
% 

figure(44)
plot(x, pdf_filtered, 'ro-')



% Build and Plot the projected RBF on the old original RBF centers
%
pdf_filtered_on_org_grid = zeros(1,N);
for i = 1:N_c


  pdf_filtered_on_org_grid = pdf_filtered_on_org_grid + new_lambdas(i).*exp(-eps_c(i)^2.*(x - x_c(i)).^2 );

end

figure(44)
hold on
plot(x, pdf_filtered_on_org_grid, 'kx-')
legend('PDF from FilterStep', 'PDF reprojected on to the original grid')

Integral_of_pdf_filtered_on_org_grid = trapz(x, pdf_filtered_on_org_grid)

%plot(x, pdf_approx, 'co-')





