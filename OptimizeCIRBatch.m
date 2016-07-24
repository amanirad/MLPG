%
%
% Optimize Likelihood estimate
%



close all
clear all
thetaSave=zeros(100,3);
timeSave=zeros(100,1);

parfor nsim=1:100
[OptVals, time] =  CIRBatch(nsim);  
    thetaSave(nsim,:)=OptVals;
    timeSave(nsim)=time;
end

