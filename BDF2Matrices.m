% Create structures with the pertinent matrices.
%[Ael,Aem,Abl,Abm,Bel,Bem]=CIRMatrices(xc,xe,xb,phi,ep,sigma_p,a_p,b_p)
function S=BDF2Matrices(beta0,Ael,Aem,Abl,Abm,Bel,Bem); 
  Cel = Ael - beta0*Bel;
  Cem = Aem - beta0*Bem; 
  %
  % LU-factorize the small square matrix
  %
  [Lb,Ub,Pb] = lu(Abm);
  %
  % Form the Schur complement
  %
  Sel = Ub\(Lb\(Pb*Abl));
  Sel = Cel - Cem*Sel;
  %
  % QR-factorize the Schur complement
  %
  [Q,R]=qr(Sel,0);

  S.Q = Q;
  S.R = R;
  S.Lb = Lb;
  S.Ub = Ub;
  S.Pb = Pb;
  S.Cem = Cem;
  S.Abl = Abl;
  S.Ael = Ael;
  S.Aem = Aem;