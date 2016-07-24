clc
a_p=2;
b_p=4;
sigma_p=1;
T=1;
x0=1;
start=0.1;
ep=1;
M=40;
N=100;
Ns=200;
rg=[0 10];
nodetype='uni';
method={'dir'};
[maxnrm,v,xs,error,xls,res,time,timerror,timeres,ao,mem]=CIRMain(a_p,b_p,sigma_p,T,x0,start,ep,M,N,Ns,rg,nodetype,method)