function XT=pTrnScl(Dsc,X,Prm);
% function XT=pTrnScl(Dsc,X,Prm);
%
% Transform between different marginal scales
% Any inf or -inf in output is replaced by NaN
%
% Possibilities are M (measured, GnrPrt), F (Frechet), L (Laplace), N (normal / Gaussian).
%
% Note that on uniform scale XUValue<0 indicates that a TAIL probability is given, which needs to be
% interpreted as GivenXU=-(1-XU) in subsequent transformations, so that XU=1+GivenXU. 
% This allows very small tail probabilities to be handled more appropriately.
%
% Dsc  DeSCription of transformation required
% X    Sample of data to be transformed
% Prm  Generalised Pareto PaRaMeters to use for the tranformation
%      This is either 2 x 1 vector (xi, sgm, with psi=0 assumed) or 3 x 1
%
% P. Jonathan 2022

%% Test decks if no input
if nargin==0;
    %
    %Check basic transformations starting on uniform scale
    tD={'F2U','U2F';'L2U','U2L';'M2U','U2M';'N2U','U2N';};
    tU1=[1e-6 0.1:0.1:0.9 1-1e-6]'*ones(1,2);
    for iD=1:4;
        if strcmp(tD{iD,1}(1),'M');
            tU2=pTrnScl(tD{iD,1},pTrnScl(tD{iD,2},tU1,[-0.2 1 3]'),[-0.2 1 3]');
        else;
            tU2=pTrnScl(tD{iD,1},pTrnScl(tD{iD,2},tU1));
        end;
        fprintf(1,'%g %g\n',mean((tU2-tU1).*2));
    end;
    %
    %Check compound tranformations starting on measured scale
    tD={'F2M','M2F';'L2M','M2L';'N2M','M2N';};
    Xi=-0.2;
    Sgm=1;
    Psi=2;
    tU1=(1:100)';
    tU1=tU1(tU1>Psi); %Admit threshold EXCEEDANCES only
    if Xi<0;
        tU1=tU1(tU1<Psi-Sgm/Xi); %Admit only points below UEP
    end;
    for iD=1:3;
        tU2=pTrnScl(tD{iD,1},pTrnScl(tD{iD,2},tU1,[Xi Sgm Psi]'),[Xi Sgm Psi]');
        fprintf(1,'%g\n',mean((tU2-tU1).*2));
    end;
    %
    return;
end;

if nargin==3; %Insert a zero as the GP threshold if not specified
    if size(Prm,1)==2; 
        Prm(3)=0;
    end;
end;

if strcmp(Dsc,'F2U'); %Frechet to Uniform
    XT=F2U(X);
    Exc=XT<0; %Tails present (so -XTOutput is actually (1-XT), and XT=1+XTOutput)
    XT(Exc)=1; %Set to unity since XTOutput is close to zero
end;

if strcmp(Dsc,'F2M'); %Frechet to Uniform to Measured (GP)
    XT=U2M(F2U(X),Prm);
end;

if strcmp(Dsc,'L2M'); %Laplace to Uniform to Measured (GP)
    XT=U2M(L2U(X),Prm);
end;

if strcmp(Dsc,'L2U'); %Laplace to Uniform
    XT=L2U(X);
    %Handle negative outputs indicating tail function used
    Exc=XT<0; %Tails present (so -XTOutput is actually (1-XT), and XT=1+XTOutput)
    XT(Exc)=1; %Set to unity since XTOutput is close to zero
end;

if strcmp(Dsc,'M2F'); %Measured to Uniform to Frechet
    XT=U2F(M2U(X,Prm));
end;

if strcmp(Dsc,'M2L'); %Measured to Uniform to Laplace
    XT=U2L(M2U(X,Prm));
end;

if strcmp(Dsc,'M2N'); %Meaasured to Uniform to Normal
    XT=U2N(M2U(X,Prm));
end;

if strcmp(Dsc,'M2U'); %Meaasured to Uniform
    XT=M2U(X,Prm);
    %Handle negative outputs indicating tail function used
    Exc=XT<0; %Tails present (so -XTOutput is actually (1-XT), and XT=1+XTOutput)
    XT(Exc)=1; %Set to unity since XTOutput is close to zero
end;

if strcmp(Dsc,'N2M'); %Normal to Uniform to Measured
    XT=U2M(N2U(X),Prm);
end;

if strcmp(Dsc,'N2U'); %Normal to Uniform
    XT=N2U(X);
    Exc=XT<0; %Tails present (so -XTOutput is actually (1-XT), and XT=1+XTOutput)
    XT(Exc)=1; %Set to unity since XTOutput is close to zero
end;

if strcmp(Dsc,'U2F'); %Uniform to Frechet
    XT=U2F(X);
end;

if strcmp(Dsc,'U2L'); %Uniform to Laplace
    XT=U2L(X);
end;

if strcmp(Dsc,'U2M'); %Uniform to Measured 
    XT=U2M(X,Prm);
end;

if strcmp(Dsc,'U2N'); %Uniform to Normal 
    XT=U2N(X);
end;

%% Ensure that any inf or -inf is replaced by NaN in output
XT(isinf(XT))=NaN;

%% Complete
return;

%% Individual functions
function XU=F2U(XF);
%F_F(x) = U = exp(-1/x)
[n,p]=size(XF);
XU=nan(n,p);
for j=1:p;
    XU(:,j)=exp(-1./XF(:,j));
    %Handle exceptions when XF is huge, resulting in MATLAB returning XU=1
    Exc=XU(:,j)==1; % XU=1 cannot be actually true!
    XU(Exc,j)=-1./XF(Exc,j);% Use negative value to indicate TAIL probability
end;
return;

function XU=L2U(XL);
[n,p]=size(XL);
XU=nan(n,p);
for j=1:p;
    XU(XL(:,j)<=0,j)=0.5*exp(XL(XL(:,j)<=0,j));
    XU(XL(:,j)>0,j)=1-0.5*exp(-XL(XL(:,j)>0,j));
    %Handle exceptions when XF is huge, resulting in MATLAB returning XU=1
    Exc=XU(:,j)==1; % XU=1 cannot be actually true!
    XU(Exc,j)=-0.5*exp(-XL(Exc,j));% Use negative value to indicate TAIL probability
end;
return;

function XU=M2U(XM,Prm);
%Generalised Pareto to Uniform
Xi=Prm(1);
Sgm=Prm(2);
Psi=Prm(3);
[n,p]=size(XM);
XU=nan(n,p);
for j=1:p;
    XU(:,j)=gpcdf(XM(:,j),Xi,Sgm,Psi);
    %Handle exceptions when XM is really close to the upper end point of the GP, indicated by Exc statement here
    Exc=XU(:,j)==1 & (1+Xi*(XM(:,j)-Psi)/Sgm).^(-1/Xi)>0; % XU=1, but handle finite upper end point
    XU(Exc,j)=-(1+Xi*(XM(Exc,j)-Psi)/Sgm).^(-1/Xi);% Use negative value to indicate TAIL probability
end;
return;

function XU=N2U(XN);
%F_N(x) = U Gaussian
[n,p]=size(XN);
XU=nan(n,p);
for j=1:p;
    XU(:,j)=normcdf(XN(:,j));
    %Handle exceptions when XN is really large, indicated by Exc statement here
    Exc=XU(:,j)==1 & normcdf(XN(:,j),'upper')>0; %Only correct when the tail prob is estimated >0
    XU(Exc,j)=-normcdf(XN(Exc,j),'upper');% Use negative value to indicate TAIL probability
end;
return;

function XF=U2F(XU);
%F_F(x) = U = exp(-1/x)
[n,p]=size(XU);
XF=nan(n,p);
for j=1:p;
    XF(:,j)=-1./log(XU(:,j));
    %Handle exceptions when XU<0 indicates TAIL probability given
    %Use MATLAB log1p function log1p(x) == log(1+x)
    Exc=XU(:,j)<0; %Tails present (so -XUInput is actually (1-XU), and XU=1+XUInput)
    XF(Exc,j)=1./log1p(-XU(Exc,j));
end;
return;

function XL=U2L(XU);
% Uniform to standard Laplace
[n,p]=size(XU);
XL=nan(n,p);
for j=1:p;
    XL(:,j)=sign(0.5-XU(:,j)).*log(2*min(1-XU(:,j),XU(:,j)));
    %Handle exceptions when XU<0 indicates TAIL probability given
    Exc=XU(:,j)<0; %Tails present (so -XUInput is actually (1-XU), and XU=1+XUInput)
    XL(Exc,j)=sign(0.5-1-XU(Exc,j)).*log(2*min(-XU(Exc,j),1+XU(Exc,j)));
end;
return;

function XM=U2M(XU,Prm);
%F_GP(x) is generalised Pareto
%Calculate on log scale for better stability
Xi=Prm(1);
Sgm=Prm(2);
Psi=Prm(3);
[n,p]=size(XU);
logXM=nan(n,p);
for j=1:p;
    if Xi<0;
        logXM(:,j)=log(Sgm)-log(Xi)+log((1-XU(:,j)).^(-Xi)-1);
    elseif Xi==0;
        logXM(:,j)=log(-log(1-XU(:,j)))+log(Sgm);
    else;
        logXM(:,j)=-Xi*log(1-XU(:,j))+log(Sgm)-log(Xi)+log(1-(1-XU(:,j)).^Xi);
    end;
    %Handle exceptions when XU<0 indicates TAIL probability given
    Exc=XU(:,j)<0; %Tails present (so -XUInput is actually (1-XU), and XU=1+XUInput)
    if Xi<0;
        logXM(Exc,j)=log(Sgm)-log(Xi)+log((-XU(Exc,j)).^(-Xi)-1);
    elseif Xi==0;
        logXM(Exc,j)=log(-log(-XU(Exc,j)))+log(Sgm);
    else;
        logXM(Exc,j)=-Xi*log(-XU(Exc,j))+log(Sgm)-log(Xi)+log(1-(-XU(Exc,j)).^Xi);
    end;
end;
%Exponentiate
XM=exp(logXM)+Psi;
return;

function XN=U2N(XU);
% XN is inverse cdf evaluated at XU;
[n,p]=size(XU);
XN=nan(n,p);
for j=1:p;
    XN(:,j)=norminv(XU(:,j));
    %Handle exceptions when XU<0 indicates TAIL probability given
    Exc=XU(:,j)<0; %Tails present
    XN(Exc,j)=sqrt(2).*erfcinv(2*(-XU(Exc,j)));
end;
return;