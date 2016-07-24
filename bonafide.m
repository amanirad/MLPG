%
% Creates a bona fide density by successive projections
% and reintegrations. A bona fide density is a density that
% is both positive definit and integrates to one.
%
% (1) 
% http://projecteuclid.org/DPubS/Repository/1.0/Disseminate?view=body&id=pdf_1&handle=euclid.aos/1176350182
%
% intput x-grid and density
function [res,error] = bonafide(y,f)

k=0;
% Algorithm on p. 1614 in (1)
fk = f;
while k<1e2;
    fk1 = max(fk, 0);
    Ck1 =trapz(y, fk1);
    if(abs(Ck1-1)<1e-10)
        break;
    end    
    % f_k+2 = f_k+1 - (C_k+1 -1)/(   h(x)* int( 1/h(x) dx) ) where h(x) is
    % a weight function. Special case : h =1 gives the scheme below=>
    fk = fk1 - ( Ck1-1)./(max(y)-min(y));
    k = k+1;
end

if(k==1e2)
error = 1;
else 
    error =0;
end

res=fk1;