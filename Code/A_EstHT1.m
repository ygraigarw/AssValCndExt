function Mdl=A_EstHT1(Dat,HTThrU,nBS);
%function Mdl=A_EstHT1(Dat,HTThrU,nBS);
%
% Associated values: estimate marginal and conditional extremes models, latter using fmincon
%
% P. Jonathan, R. Towe 2022

SmpSiz=size(Dat.M,1);

%Loop over bootstrap resamples
for iBS=1:nBS; 
        
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
    tGPPrmX=gpfit(YBS(:,1))';
    tGPPrmY=gpfit(YBS(:,2))';
    %Add to structure
    Mdl.GP.Prmh(:,1,iBS)=tGPPrmX;
    Mdl.GP.Prmh(:,2,iBS)=tGPPrmY;
    
    %% Transform to standard Laplace
    tL1=pTrnScl('M2L',YBS(:,1),[tGPPrmX;0]);
    tL2=pTrnScl('M2L',YBS(:,2),[tGPPrmY;0]);
    
    %% Estimate HT threshold on Laplace scale
    HTThrL=pTrnScl('U2L',HTThrU);
    
    %% Fit Heffernan and Tawn - estimate parameters + residuals
    Prm0=A_HTStrtSln(tL1,tL2,0);
    A=[];b=[];Aeq=[];beq=[];lb=[1e-6;-100;-100;1e-6];ub=[1.2-1e-6;1-1e-6;100;100];
    nonlcon=[];options=optimoptions('fmincon','Display', 'off');
    Prmh=fmincon(@(Prm)A_HTMLE(Prm,[tL1 tL2],HTThrL),Prm0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    
    %% Estimate residuals
    x=tL1;
    Kep=x>HTThrL;
    x=tL1(Kep);
    y=tL2(Kep);
    Rsd=(y-Prmh(1)*x-Prmh(3)*x.^Prmh(2))./(Prmh(4)*x.^Prmh(2));
    %Add to structure
    Mdl.HT.Prmh(:,iBS)=Prmh;
    Mdl.HT.Rsd{iBS}=Rsd;
    
end;

%% Complete
return;