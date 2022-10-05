function [CndExp,CndMdn]=A_EYMgXF(XFCnd,Alp,PrmGP);
% function [CndExp,CndMdn]=A_EYMgXF(XFCnd,Alp,PrmGP);
%
% Conditional expectation on MEASURED SCALE of the bivariate logistic (also called bivariate Gumbel)
% Conditioning value XFCnd on Frechet scale
%
% P. Jonathan, R. Towe 2022

%% Test case
if nargin==0;
    %XFCnd=100;
    XFCnd=10000;%1e20;%-1/log(0.999999999999);
    Alp=0.1;%0.95;
    PrmGP=[0.2;1];
end;

%% Set-up
% Set-up parameters (probably don't need to change)
% These parameters should be OK, checked for Alp in range 0.01 to 0.99 and XFCnd<=100

%% Grid locations
nDlt=1000;
Lwr=-5;
Upp=log(XFCnd*10000);
LgrYF=linspace(Lwr,Upp,nDlt);
DltLgrF=diff(LgrYF(1:2));
YF=exp(LgrYF);

%% Estimate conditional density for bivariate logistic
CndDns=A_fYFgXF(YF,XFCnd,Alp);

%% Estimate conditional expectation
Wgh=pTrnScl('F2M',YF,PrmGP);
tScl=sum(CndDns.*YF)*DltLgrF; %scale in case density does not sum to unity
CndDns=CndDns/tScl; %correct the conditional density to sum to unity
CndExp=sum(Wgh.*CndDns.*YF)*DltLgrF;

%% Estimate conditional median
tCdf=cumsum(CndDns.*YF)*DltLgrF/tScl;
tLct=find(tCdf>=0.5);
if isempty(tLct)==0;
    MdnYFgXF=YF(tLct(1));
    CndMdn=pTrnScl('F2M',MdnYFgXF,PrmGP);
else;
    CndMdn=NaN;
end;

if nargin==0;
    clf;
    fprintf(1,'Integral of uncorrected CndDns (using log scale) = %g\n',tScl);
    fprintf(1,'Integral of corrected CndDns (using log scale) = %g\n',sum(CndDns.*YF)*DltLgrF);
    fprintf(1,'Conditioning value (measurement scale) = %g\n',pTrnScl('F2M',XFCnd,PrmGP));
    fprintf(1,'CndExp (using log scale, on measurement scale) = %g\n',CndExp);
    fprintf(1,'CndMdn (using log scale, on measurement scale) = %g\n',CndMdn);
    subplot(2,2,1);plot(log10(YF),CndDns,'k-');
    subplot(2,2,2);plot(log10(YF),Wgh,'k-');
    subplot(2,2,3);plot(log10(YF),Wgh.*CndDns,'k-');
    subplot(2,2,4); hold on;
    plot(log10(YF),cumsum(CndDns.*YF)*DltLgrF,'k-');
end;

%% Complete
return;