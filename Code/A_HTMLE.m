function nll=A_HTMLE(Prm,DatL,ThrL);

% Isolate threshold exceedances
x=DatL(:,1);
Kep=x>ThrL;
x=DatL(Kep,1);
y=DatL(Kep,2);

% Parameter values
a=Prm(1);
b=Prm(2);
mu=Prm(3);
sgm=Prm(4);

% Negative log Gaussian likelihoood
x2b=x.^b;
sx2b=sgm*x.^b;
%nll=sum(log(sgm*x.^b)+0.5*((y-a*x-mu*x.^b)./(sgm*x.^b)).^2);
nll=sum(log(sx2b)+0.5*((y-a*x-mu*x2b)./(sx2b)).^2);

return;