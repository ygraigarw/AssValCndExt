function Mdl=A_EstCplTruMrg(Dat,Dpn,nBS,Xix,Xiy);
%function Mdl=A_EstCpl(Dat,Dpn,nBS);
%
% Associated values: estimate copula models with known margins
%
% P. Jonathan, R. Towe 2022

SmpSiz=size(Dat.M,1);

if strcmp(Dpn,'Gss')==1; %Gaussian dependence
    
    for iBS=1:nBS; %Loop over bootstraps
        
        %% Bootstrap indicator
        if iBS==1;
            BSInd=(1:SmpSiz)';
        else;
            BSInd=randsample((1:SmpSiz)',SmpSiz,1);
        end;
        %Add to structure
        Mdl.BSInd(:,iBS)=BSInd;
        
        %% Bootstrap resample
        YBS=Dat.M(BSInd,:);
        
        %% Marginal models
        tGPPrmX=[Xix 1]'; %true margin
        tGPPrmY=[Xiy 1]'; %true margin
        %Add to structure
        Mdl.GP.Prmh(:,1,iBS)=tGPPrmX;
        Mdl.GP.Prmh(:,2,iBS)=tGPPrmY;
        
        %% Transform to standard Normal
        tN1=pTrnScl('M2N',YBS(:,1),[tGPPrmX;0]);
        tN2=pTrnScl('M2N',YBS(:,2),[tGPPrmY;0]);
        
        %% Fit correct copula dependence to full sample
        Prm0=0.5;
        A=[];b=[];Aeq=[];beq=[];lb=[0];ub=[1];
        nonlcon=[];options=optimoptions('fmincon','Display','off');
        Prmh=fmincon(@(Prm)BvrGssNLL(Prm,[tN1 tN2]),Prm0,A,b,Aeq,beq,lb,ub,nonlcon,options);
        %Add to structure
        Mdl.Cpl.Prmh(:,iBS)=Prmh;
        
    end;
    
elseif strcmp(Dpn,'Lgs')==1; %Logistic dependence
    
    for iBS=1:nBS; %Loop over bootstraps
        
        %% Bootstrap indicator
        if iBS==1;
            BSInd=(1:SmpSiz)';
        else;
            BSInd=randsample((1:SmpSiz)',SmpSiz,1);
        end;
        %Add to structure
        Mdl.BSInd(:,iBS)=BSInd;
        
        %% Bootstrap resample
        YBS=Dat.M(BSInd,:);
        
        %% Marginal models
        tGPPrmX=[Xix 1]'; %true margin
        tGPPrmY=[Xiy 1]'; %true margin
        %Add to structure
        Mdl.GP.Prmh(:,1,iBS)=tGPPrmX;
        Mdl.GP.Prmh(:,2,iBS)=tGPPrmY;
        
        %% Transform to standard Frechet
        tF1=pTrnScl('M2F',YBS(:,1),[tGPPrmX;0]);
        tF2=pTrnScl('M2F',YBS(:,2),[tGPPrmY;0]);
        
        %% Fit correct copula dependence to full sample
        Prm0=0.5;
        A=[];b=[];Aeq=[];beq=[];lb=[0];ub=[1];
        nonlcon=[];options=optimoptions('fmincon','Display','off');
        Prmh=fmincon(@(Prm)BvrLgsNLL(Prm,[tF1 tF2]),Prm0,A,b,Aeq,beq,lb,ub,nonlcon,options);
        %Add to structure
        Mdl.Cpl.Prmh(:,iBS)=Prmh;
        
    end;
    
end;

%% Complete
return;

%% Bivariate Gaussian negative log likelihood
function NLL=BvrGssNLL(Rho,X);

n=size(X,1);
InvSgm=(1/(1-Rho^2))*[1 -Rho;-Rho 1];

%NLL=n*log(2*pi)+(n/2)*log(1-Rho^2)+0.5*sum(diag(X*InvSgm*X'));
NLL=n*log(2*pi)+(n/2)*log(1-Rho^2)+0.5*sum(X.*(InvSgm*X')',[1 2]);

%% Complete
return;

%% Bivariate logistic negative log likelihood
function NLL=BvrLgsNLL(Alp,X);

% G=(X(:,1).^(-1/Alp)+X(:,2).^(-1/Alp));
% V=G.^Alp;
% 
% dGdx=(-1/Alp)*X(:,1).^(-1/Alp-1);
% dGdy=(-1/Alp)*X(:,2).^(-1/Alp-1);
% 
% dVdx=Alp*G.^(Alp-1).*dGdx;
% dVdy=Alp*G.^(Alp-1).*dGdy;
% 
% d2Vdxdy=Alp*(Alp-1)*G.^(Alp-2).*dGdx.*dGdy;
% 
% NLL=sum(V-log(dVdx.*dVdy-d2Vdxdy));

% G=(X(:,1).^(-1/Alp)+X(:,2).^(-1/Alp));
% GAM2=G.^(Alp-2);
% GAM1=GAM2.*G;
% V=GAM1.*G;
% 
% dGdx=(-1/Alp)*X(:,1).^(-1/Alp-1);
% dGdy=(-1/Alp)*X(:,2).^(-1/Alp-1);
% 
% dVdx=Alp*GAM1.*dGdx;
% dVdy=Alp*GAM1.*dGdy;
% 
% d2Vdxdy=Alp*(Alp-1)*GAM2.*dGdx.*dGdy;
% 
% NLL=sum(V-log(dVdx.*dVdy-d2Vdxdy));

XPAM1=X.^(-1/Alp-1);
%G=(X(:,1).^(-1/Alp)+X(:,2).^(-1/Alp));
G=sum(XPAM1.*X,2);
GAM2=G.^(Alp-2);
GAM1=GAM2.*G;
V=GAM1.*G;

dGdx=(-1/Alp)*XPAM1(:,1);
dGdy=(-1/Alp)*XPAM1(:,2);

dVdx=Alp*GAM1.*dGdx;
dVdy=Alp*GAM1.*dGdy;

d2Vdxdy=Alp*(Alp-1)*GAM2.*dGdx.*dGdy;

NLL=sum(V-log(dVdx.*dVdy-d2Vdxdy));

%% Complete
return;