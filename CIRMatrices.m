function [Ael,Aem,Abl,Abm,Bel,Bem]=CIRMatrices(xc,xb,xe,phi,ep,sigma_p,a_p,b_p)
%
% This function returns the operator matrix for the CIR-PDE in one
% dimension. This includes the boundary condition equations as well.
% As the interpolation matrix is computed in the process, we return that
% as well. Maybe we can't separate the operator and the rest...  
% Input xe are assumed to be interior points only or all the points where
% the least squares are applied. Could be also the boundary points
% actually. But if both are included, the matrix is closer to singular in
% the square case.  
% The boundary points are used for the "boundary columns"  
  
  Ne = size(xe,1);
  N = size(xc,1);
  Nb = size(xb,1);
  %
  % Form some coefficients for the operator
  %
  sig2=sigma_p*sigma_p';
  X=spdiags(xe,0,Ne,Ne);
  I=speye(Ne);
  %
  % Distances between all types of nodes
  %
  r_ec = xcdist(xe,xc,1);
  r_bc = xcdist(xb,xc,1);
  r_bb = xcdist(xb,xb,1);
  r_eb = xcdist(xe,xb,1);
  %
  % RBF interpolation matrices
  %
  Ael = RBFmat(phi,ep,r_ec,'0');
  Aem = RBFmat(phi,ep,r_eb,'0');
  Abl = RBFmat(phi,ep,r_bc,'0');
  Abm = RBFmat(phi,ep,r_bb,'0');
  %
  % First derivatives
  %
  dim = 1;
  A1el = RBFmat(phi,ep,r_ec,'1',dim);
  A1em = RBFmat(phi,ep,r_eb,'1',dim);
  A1bl = RBFmat(phi,ep,r_bc,'1',dim);
  A1bm = RBFmat(phi,ep,r_bb,'1',dim);
  %
  % Second derivatives
  %
  A2el = RBFmat(phi,ep,r_ec,'2',dim);
  A2em = RBFmat(phi,ep,r_eb,'2',dim);
  %
  % Form the specific CIR operator matrices 
  % Lf = 0.5*x*sig2*f'' + (a_p(x-b)+sig2)*f' + a_p*f
  %
  Bel = 0.5*sig2*X*A2el + (a_p*(X-b_p*I)+sig2*I)*A1el + a_p*Ael;
  Bem = 0.5*sig2*X*A2em + (a_p*(X-b_p*I)+sig2*I)*A1em + a_p*Aem;
  %
  % Flux conserving boundary conditions
  % The first boundary is Dirichlet, the second Robin: 
  % (a(x-b)+sig2/2)*f+sig2/2*x*f'
  %
  Abl(2,:) = (a_p*( xb(2) - b_p) + sig2/2 )*Abl(2,:) ...
      + sig2/2*xb(2)*A1bl(2,:);
  Abm(2,:) = (a_p*( xb(2) - b_p) + sig2/2 )*Abm(2,:) ...
      + sig2/2*xb(2)*A1bm(2,:);    

