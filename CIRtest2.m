
phi='gs';

a=1;
b=2;
s0=1;
sig0=0.5;
mu = @(t) a*(b-s0)*t;
sigma = @(t) sig0*sqrt(s0*t);

t1=1e-3;
T=0.5;
%
% Shape parameter values for the first and the last time
%
ep(1)=1/sqrt(2)/sigma(t1);
ep(2)=1/sqrt(2)/sigma(T);

%ep=ep/4;
%
% How many nodes should we use to have phi(h/2)=1/2?
%
h=2*sqrt(log(1/0.7))./ep;

L=4;
rg=[0 L];
N=10+round(L./h+1);

% Time steps.
M= 500;

% Evaluation points
Ns=100;

ctype='uni';
etype='uni';
method={'dir'};

show='yes';
col='k';
%
% With 0.7 it goes down a bit. Saturation or conditioning?
s=0.9;

Nvec=[100:16:200];
Nvec=100;
% 
% Computethe corresponding epsilon
% 
d = L./(Nvec-1);

epvec=2*sqrt(-log(s))./d;

for k=1:length(Nvec)
  N=Nvec(k);
  ep=epvec(k);
%
% Least squares points.
%
Nls=4*max(N);
% Rung CIRMLMain instead of LeastSquares.
maxnrm(k)=CIRcompare(a,b,sig0,T,s0,t1,ep,M,N,Ns,rg,ctype,method,Nls);
end
%
% Should we keep the integral = 1 or the integral over the domain we have
% to constant?
%
% QUESTION: How do we do the right hand side with the integral if we have
% multiple levels?
% Can our current code handle multiple levels?
% Somehtin to think through.
% We try to get the integral to one and then we get a residual, but we
% don't get a residual because boundary conditions are satisfied exactly.

%
% L?sningen driver ?t fel h?ll eller inte alls.
%