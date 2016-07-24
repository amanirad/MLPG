SetPaths
phi='gs';
M=500;
N=200;
Ns=200;
rg=[0 4];
ctype='uni';
etype='uni';
method={'dir'}
show='no';
epvec=logspace(0,1,50);
for k=1:length(epvec)
  ep=epvec(k);
  [maxnrm(k),finnrm,xs,error,xls,res,time,timerror,timeres,ao,mem]=BSMain1D(phi,ep,M,N,Ns,rg,ctype,etype,method,show,'b');
end
pos = find(maxnrm<1);
semilogy(epvec(pos),maxnrm(pos),'b')

epres = epvec(pos);
nrmres=maxnrm(pos);