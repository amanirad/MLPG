function [logp,logp2]=CIRparticleFilter(data,t,param,K,M)
% 
% data          row vector
% t             row vector
% param = [a b sig]
% K # nbr of particles
% M # Nbr of subsets
% Gamma = 1e-5  Measurement error variance

N=length(t);
Dt=t(2)-t(1);
dt=Dt/M;
a=param(1);
b=param(2);
sig=param(3);
MyGamma=1e-5;
logp=zeros(1,N);

% Set up initial conditions
xif=data(1);
wf=ones(K,1)/K;

for n=1:N
    % predict
    xi=xif;
    dW=sqrt(dt)*randn(K,M-1);
    for m=1:(M-1)
        xi=abs(xi+a*(b-xi)*dt+sig*sqrt(xi).*dW(:,m));
    end
    % Better computation of the log-likelihood
    logp(n)=log(sum(wf.*normpdf(data(n),xi+a*(b-xi)*dt,sqrt(sig^2*xi*dt+MyGamma))));

    xipr=xi+a*(b-xi)*dt+sig*sqrt(dt)*sqrt(xi).*randn(K,1);
    wpr=wf.*normpdf(data(n),xipr,sqrt(MyGamma));
    logp2(n)=log(sum(wpr));
    
    sw=sum(wpr);
    if sw>0
        w=wpr/sw;
    else
        w=ones(K,1)/K;
    end
    
    % resample
    I=randsample(1:K,K,'true',w);
    xif=xipr(I);
end
logp=sum(logp);
logp2=sum(logp2);