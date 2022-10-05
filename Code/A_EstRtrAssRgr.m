function VlsEst=A_EstRtrAssRgr(Mdl,AnnRat,RtrPrd);
%function VlsEst=A_EstRtrAssRgr(Mdl,AnnRat,RtrPrd);
%
% Associated values: estimate return and associated values using copula models
%
% P. Jonathan, R. Towe 2022

VlsEst.Rtr=nan(1,5); %Marginal return values
VlsEst.Ass=nan(5,4); %Associated values

nBS=size(Mdl.GP.Prmh,3); %Number of bootstrap resamples

Rtr.U=1-(-1/AnnRat)*log(1-1/RtrPrd); %know this in closed form
Rtr.L=pTrnScl('U2L',Rtr.U);

%% Return values

%q1
ExpPrmh=mean(permute(Mdl.GP.Prmh(:,1,:),[3 1 2]))';
VlsEst.Rtr(1)=pTrnScl('U2M',Rtr.U,[ExpPrmh(1);ExpPrmh(2);0]); %transform to physical (measured) scale

% q2
tQnt=nan(nBS,1);
for iBS=1:nBS;
    tQnt(iBS)=pTrnScl('U2M',Rtr.U,[Mdl.GP.Prmh(1,1,iBS);Mdl.GP.Prmh(2,1,iBS);0]);
end;
VlsEst.Rtr(2)=mean(tQnt);

% q3
DInc=100;DMxm=10000;SlnLwr=0;
for iI=1:4; %loop to refine estimate
    D=(0:DInc:DMxm-DInc)'+SlnLwr; nD=size(D,1);
    tCdf=A_AnnCdf(D,permute(Mdl.GP.Prmh(1,1,:),[3 2 1]),permute(Mdl.GP.Prmh(2,1,:),[3 2 1]),AnnRat,1)';
    SlnLwr=D(tCdf<1-1/RtrPrd);SlnLwr=SlnLwr(end);
    DInc=DInc/100; DMxm=DMxm/100;
end;
VlsEst.Rtr(3)=SlnLwr;

% q4
DInc=100;DMxm=10000;SlnLwr=0;
for iI=1:4; %loop to refine estimate
    D=(0:DInc:DMxm-DInc)'+SlnLwr; nD=size(D,1);
    tCdf=A_AnnCdf(D,permute(Mdl.GP.Prmh(1,1,:),[3 2 1]),permute(Mdl.GP.Prmh(2,1,:),[3 2 1]),AnnRat,RtrPrd)';
    SlnLwr=D(tCdf<exp(-1));SlnLwr=SlnLwr(end);
    DInc=DInc/100; DMxm=DMxm/100;
end;
VlsEst.Rtr(4)=SlnLwr;

% estimate q5
% calculate quantile of annual max distribution for each bootstrap
% take median over estimated quantiles
tQnt=nan(nBS,1);
for iBS=1:nBS;
    tQnt(iBS)=pTrnScl('U2M',Rtr.U,[Mdl.GP.Prmh(1,1,iBS);Mdl.GP.Prmh(2,1,iBS);0]);
end;
VlsEst.Rtr(5)=median(tQnt);

%% Associated values on M scale using definitions 1 and 2
for k=1:5; %Loop over return value for X
    
    %Conditioning value measured-scale
    CndValM=VlsEst.Rtr(k);
        
    EYM=nan(nBS,1);
        
    for iBS=1:nBS; %Loop over bootstraps
        
        %Linear regression estimate of EY | X = RV on measurement scale
        EYM(iBS)=Mdl.Cpl.Prmh(1,1,iBS) + CndValM*Mdl.Cpl.Prmh(2,1,iBS);
        
    end; % bootstrap
    
    %Manage Inf
    EYM(isinf(EYM))=NaN;
        
    %We are interested in mean and median over bootstraps
    VlsEst.Ass(k,1)=mean(EYM(~isnan(EYM)));
    VlsEst.Ass(k,2)=median(EYM(~isnan(EYM)));
    
    %Record problematic bootstraps
    VlsEst.AssOK(k)=sum(~isnan(EYM));
    
end; % choice of return value

%% q_23, q_24
for iBS=1:nBS;
    
    %Return value for X on U, L and M scales
    Rtr.U=1-(-1/AnnRat)*log(1-1/RtrPrd); %know this in closed form
    Rtr.L=pTrnScl('U2L',Rtr.U);
    Rtr.M(iBS)=pTrnScl('U2M',Rtr.U,[Mdl.GP.Prmh(:,1,iBS);0]);
    
    %Linear regression estimate of EY | X = RV on measurement scale
    EYM(iBS)=Mdl.Cpl.Prmh(1,1,iBS) + Rtr.M(iBS)*Mdl.Cpl.Prmh(2,1,iBS);
    
end;

%Manage Inf
EYM(isinf(EYM))=NaN;

%We are interested in means and medians over bootstraps
VlsEst.Ass(2,3)=mean(EYM(~isnan(EYM)));
VlsEst.Ass(2,4)=median(EYM(~isnan(EYM)));
VlsEst.Ass(5,3)=VlsEst.Ass(2,3);
VlsEst.Ass(5,4)=VlsEst.Ass(2,4);

%Record any problematic bootstraps
VlsEst.AssOK(6)=sum(~isnan(EYM));

%% Complete
return;