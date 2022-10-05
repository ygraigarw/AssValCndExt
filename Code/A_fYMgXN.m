function CndDns=A_fYMgXN(YN,XNCnd,Rho,PrmGP);
%function CndDns=fYMgXN(YN,XNCnd,Rho,PrmGP);
%
%Calculate transformed conditional density for Gaussian.
%Need this for the g_Y^{-1}(Y_L|X_L) calculation.

%% Grid specification
DltN=0.01;
YNMxm=100;

if nargin==0; % Test
    YN=(-YNMxm:DltN:YNMxm)';
    XNCnd=5;
    Rho=0.1;
    PrmGP=[0;1];
end;

%% Parameters for conditional Gaussian
n=size(YN,1);
Mu=Rho*XNCnd;
Tau=sqrt(1-Rho^2); % NB use Tau for standard deviation, reserve Sgm for GP scale

%% Check that function tails off appropriately on log scale
t1=-YNMxm:DltN:YNMxm;
tU=normcdf(t1); 
%Transform to measured scale
Xi=PrmGP(1);
Sgm=PrmGP(2);
t2a=log(pTrnScl('U2M',tU,PrmGP));
t2b=pLgrNrmDns(t1,[Mu;Tau]);
t2=t2a+t2b;
yInt=find(t2>-100 & isinf(t2)==0); % this is the largest value of YN for which output>1e-100

%% Check that yInt choice looks sensible
tMxm=max(t2(yInt));
if yInt(1)==1;
    tEndMnm=t2(yInt(end));
else;
    tEndMnm=max(t2([yInt(1) yInt(end)]));
end;
Rat=tMxm-tEndMnm;
if Rat<10;
    fprintf(1,'ERROR IN CndDns. Rat=%g\n',Rat);
    CndDns=nan(n,1);
    clf; hold on;
    plot(t1,t2,'ko-');  
    plot(t1,t2a,'ro-');  
    plot(t1,t2b,'go-');  
    return;
end;

%% Estimate transformed conditional density
CndDns=zeros(n,1);
yOK=YN>t1(yInt(1)) & YN<=t1(yInt(end));
tU=normcdf(YN(yOK)); 
if Xi<0;
    tlogM=log((Sgm/Xi)*( (1-tU).^(-Xi) - 1) + 0);
elseif Xi==0;
    tlogM=log(-log(1-tU)*Sgm+0);
else;
    tlogM=-Xi*log(1-tU)+log(Sgm/Xi)+log(1-(1-tU).^Xi);
end;
CndDns(yOK)=exp(tlogM+pLgrNrmDns(YN(yOK),[Mu;Tau]));

if nargin==0;
    clf;
    plot(YN,CndDns,'k-');
end;

%% Complete
return;
