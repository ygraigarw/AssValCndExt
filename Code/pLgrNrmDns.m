function y=pLgrNrmDns(x,Prm);
% function y=pLgrNrmDns(x,Mu,Tau);
%
% log density of normal distribution with parameters Prm

Mu=Prm(1);
Tau=Prm(2); %Tau not Sgm
y=-(1/2)*log(2*pi*Tau^2)-(1/2)*((x-Mu)./Tau).^2;

return;