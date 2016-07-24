function logpSurface=CIRparticleFitlerSurface(data,t,K,M)
b=0.05;
aVec=linspace(3,7,15);
sigVec=linspace(.05,.25,20);

logpSurface=zeros(length(aVec),length(sigVec));

for na=1:length(aVec)
    aa=aVec(na)
    
    for nsig=1:length(sigVec)
        sig=sigVec(nsig);
        
        % Define local parameter
        param=[aa b sig];
        logpSurface(na,nsig)=CIRparticleFilter(data,t,param,K,M);
    end
end
