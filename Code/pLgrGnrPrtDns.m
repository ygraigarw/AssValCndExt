function logDns=pLgrGnrPrtDns(X,Prm);
% function logDns=pLgrGnrPrtDns(X,Prm);
%
% log density of the generalised Pareto distribution 

Xi=Prm(1);
Sgm=Prm(2);

if Xi<0;
    n=size(X,1);
    IsOK=X<-Sgm/Xi;
    logDns=nan(n,1);
    logDns(IsOK)=-log(Sgm)-(1+1/Xi)*log(1+(Xi/Sgm).*X(IsOK));
elseif Xi==0;
    logDns=-log(Sgm)-X/Sgm;
else;
    logDns=-log(Sgm)-(1+1/Xi)*log(1+(Xi/Sgm).*X);    
end;

return;