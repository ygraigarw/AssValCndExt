function CndExp=A_EYMgXN2(XNCnd,Rho,PrmGP);
% function CndExp=EYMgXN2(XNCnd,Rho,PrmGP);
% 
% Calculate transformed conditional expectation for Gaussian.

if nargin==0;
    XNCnd=5;
    Rho=0.9;
    PrmGP=[0.2;1];
end;

%% Grid locations
XMMxm=1000;
DltM=0.01;
XM=(DltM:DltM:XMMxm)';

%% log GP density
logfm=pLgrGnrPrtDns(XM,PrmGP);

%% log ratio of conditioned to unconditioned Gaussians
% Uses XN not XM
XN=pTrnScl('M2N',XM,PrmGP); %transform XM to XN
logRat=-0.5*log(1-Rho^2) + (XNCnd^2)/2 -0.5*((XN-XNCnd/Rho).^2)*(Rho^2)/(1-Rho^2);
logRat(isinf(logRat))=nan;

%% log integrand
logInt=log(XM)+logRat+logfm;

%% integral
CndExp=nansum(exp(logInt))*DltM;

if nargin==0;   
    clf; 
    subplot(1,3,1); plot(XM, logfm); title('log density (measurement scale)');
    subplot(1,3,2); plot(XM, XN); title('log values (normal scale)');
    subplot(1,3,3); plot(XN, logRat); title('log ratio of normal densities');
end;

return;


