function CndExp=A_EYMgXN1(XNCnd,Rho,PrmGP);
% function CndExp=EYMgXN1(XNCnd,Rho,PrmGP);
% 
% Calculate transformed conditional expectation for Gaussian.

if nargin==0; % Test
    XNCnd=5;
    Rho=0.9;
    PrmGP=[0.2;1];
end;

%% Set-up
% Set-up parameters (probably don't need to change)
% These parameters should be OK, checked for Alp in range 0.01 to 0.99 and XFCnd<=100
XNMxm=100;
DltN=0.01;
% Grid locations
XN=-XNMxm:DltN:XNMxm;
% Parameters for conditional Gaussian
Mu=Rho*XNCnd;
Tau=sqrt(1-Rho^2); %Tau not Sgm

%% Estimate conditional density
logfXNgXN=pLgrNrmDns(XN,[Mu;Tau]);

%% Estimate conditional expectation
logWgh=log(pTrnScl('N2M',XN,PrmGP));
logInt=logfXNgXN+logWgh;
logInt(isinf(logInt)==1)=nan;
CndExp=nansum(exp(logInt))*DltN;

if nargin==0;   
    clf; 
    subplot(1,3,1); plot(XN, logWgh); title('log weights');
    subplot(1,3,2); plot(XN, logfXNgXN); title('log conditional density (normal scale)');
    subplot(1,3,3); plot(XN, logInt); title('log integrand');
end;

return;
