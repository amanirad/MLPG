% NOT: I CIRMLMatrices en parameter Bugg
%
% OBSERVATIONS:
% The ep=ep(1)*2, h med 0.7 det finns negativa coefficienter. ?verlag ?r
% coefficienterna en bra approximation av l?sningen. Det som h?nder ?r
% oscillationer i svansen med varannan negativ coefficient.
% Coefficienterna avtar fint mot kanten. Hur kan man kombinera en positiv
% och en negativ p? b?sta s?tt? Justera om mha Taylor?
% Storleken p? koefficienterna minskar med tiden.
%
% Med 0.5 i h blir approximationen d?lig. "Gropar". 
% 0.9: ser bra ut, felet v?xer med tiden till 0.27
% 0.7 n?stan exakt samma fel.
% 0.99 n?stan exakt samma fel igen. 176 pkter mot 31 punkter med samma
% ep. S? vad ?r felet? Det ?r samma varje g?ng.
% Prova att r?kna ut negativa massan, men borde vara mindre...

%
% Things to find out:
% 1. Can we run multilevel
%    Yes. In one example there was a gain with two levels compared with
%    one. When we start at t=1e-2, using ep(1) we get an error of 1e-4.
%    using N=54 points. Instead using ep(2) and N=18, we get a large
%    initial error which is not completely damped out. The oscillations
%    go away, but the negative lobe stays longer. Using both grids does
%    not lead to a smaller error, but we are already at the numerical limit.  

% 2. How does convergence look?
% 3. How do we best choose parameters?
%    The exact point where it goes bad seems to be when RBFs are further
%    apart than the half point. This makes sense. The best result is at
%    the 0.995 point. Lets do the same thing with another N. Same for
%    going bad, 0.83 for best. Conditioning involved is a reasonable
%    guess. Same minimum error in both cases.
%    Varying N goes down until 0.75 approximately, but oscillates after
%    that and perhaps becomes ill-conditioned.
%    N is an even multiple of 4 works best. Why? Because of the interval
%    length, initial location of peak?



% 4. How spikey can we be?
%    Starting at 1e-3 gives a bit larger error. And more points.
%    err=0.03 with standard setting. Making ep half and N about half too
%    lands the error at 7.8e-4. 4 times gives 1.2e-3.

% 5. Can we use variable shape and location?  
phi='gs';

%
% Try to use epsilon-values that correspond to initial and final time
%
% Formula from Josef: N(a(b-s0)dt,sigma*sqrt(s0*dt))
%
% Insert parameters that we are using. This can go wrong if code is changed.
%
a=0.4;
b=0.52;
s0=1;
sig0=0.0732;
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
M= 100;

% Evaluation points
Ns=100;

ctype='uni';
etype='uni';
method={'dir'};

show='yes';
col='k';
%
% Run a range
%
%N = 70;
%epvec=[2 5 10 15 20 25 30 35 40 45];

% Choose a sweet point
% With sweet=0.5, 0.6, the error grows
% With 0.7 it goes down a bit. Saturation or conditioning?
s=0.9;

Nvec=[200:16:300];
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
maxnrm(k)=CIRMain(a,b,sig0,T,s0,t1,ep,M,N,Ns,rg,ctype,method,Nls);
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