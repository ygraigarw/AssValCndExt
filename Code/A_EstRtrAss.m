function VlsEst=A_EstRtrAss(Mdl,VlsTru,AnnRat,RtrPrd);
%function VlsEst=A_EstRtrAss(Mdl,VlsTru,AnnRat,RtrPrd);
%
% Associated values: estimate return and associated values from fitted models
%
% P. Jonathan, R. Towe 2022

nBS=size(Mdl.GP.Prmh,3); %number of bootstrap resamples

%% Estimate marginal return values
VlsEst.Rtr=nan(1,5);

%Estimate q1
% find bootstrap mean GP parameters for X
% calculate quantile of annual max distribution
ExpPrmh=mean(permute(Mdl.GP.Prmh(:,1,:),[3 1 2]))';
VlsEst.Rtr(1)=pTrnScl('U2M',VlsTru.Rtr.U,[ExpPrmh(1);ExpPrmh(2);0]); %transform to physical (measured) scale

%Estimate q2
% calculate quantile of annual max distribution for each bootstrap
% take mean over estimated quantiles
tQnt=nan(nBS,1);
for iBS=1:nBS;
    tQnt(iBS)=pTrnScl('U2M',VlsTru.Rtr.U,[Mdl.GP.Prmh(1,1,iBS);Mdl.GP.Prmh(2,1,iBS);0]);
end;
VlsEst.Rtr(2)=mean(tQnt);

%Estimate q3
% calculate predictive distribution of annual maximum
% read off 1-1/RtrPrd quantile
DInc=100;DMxm=10000;SlnLwr=0;
for iI=1:4; %loop to refine estimate
    D=(0:DInc:DMxm-DInc)'+SlnLwr; nD=size(D,1);
    tCdf=A_AnnCdf(D,permute(Mdl.GP.Prmh(1,1,:),[3 2 1]),permute(Mdl.GP.Prmh(2,1,:),[3 2 1]),AnnRat,1)';
    SlnLwr=D(tCdf<1-1/RtrPrd);SlnLwr=SlnLwr(end);
    DInc=DInc/100; DMxm=DMxm/100;
end;
VlsEst.Rtr(3)=SlnLwr;

%Estimate q4
% calculate predictive distribution of maximum over RtrPrd years
% read off exp(-1) quantile
DInc=100;DMxm=10000;SlnLwr=0;
for iI=1:4; %loop to refine estimate
    D=(0:DInc:DMxm-DInc)'+SlnLwr; nD=size(D,1);
    tCdf=A_AnnCdf(D,permute(Mdl.GP.Prmh(1,1,:),[3 2 1]),permute(Mdl.GP.Prmh(2,1,:),[3 2 1]),AnnRat,RtrPrd)';
    SlnLwr=D(tCdf<exp(-1));SlnLwr=SlnLwr(end);
    DInc=DInc/100; DMxm=DMxm/100;
end;
VlsEst.Rtr(4)=SlnLwr;

%Estimate q5
% calculate quantile of annual max distribution for each bootstrap
% take median over estimated quantiles
tQnt=nan(nBS,1);
for iBS=1:nBS;
    tQnt(iBS)=pTrnScl('U2M',VlsTru.Rtr.U,[Mdl.GP.Prmh(1,1,iBS);Mdl.GP.Prmh(2,1,iBS);0]);
end;
VlsEst.Rtr(5)=median(tQnt);

%% Estimate associated values under fitted model (see paper for definitions)
VlsEst.Ass=nan(5,4); % 5 X-return values, max of 4 associated values

%% First estimate associated values 1-2, where the X return value is precomputed
for k=1:5; %Loop over return value for X
    
    %Conditioning value measured-scale
    CndValM=VlsEst.Rtr(k);
    
    EYM=nan(nBS,1);
    
    for iBS=1:nBS;
        
        tUEP=-Mdl.GP.Prmh(2,1,iBS)/Mdl.GP.Prmh(1,1,iBS);
        if Mdl.GP.Prmh(1,1,iBS)<0;
            tOK=tUEP>CndValM;
        else;
            tOK=1;
        end;
        
        if tOK==1;
            
            %Access HT sample on Laplace scale corresponding to CndValM
            CndValL=pTrnScl('M2L',CndValM,[Mdl.GP.Prmh(:,1,iBS);0]);
            Prmh=Mdl.HT.Prmh(:,iBS);
            Rsdh=Mdl.HT.Rsd{iBS};
            %
            YL=A_YLSmp(CndValL,Prmh,Rsdh);
            %
            YM=pTrnScl('L2M',YL,[Mdl.GP.Prmh(:,2,iBS);0]);
            EYM(iBS)=mean(YM); %mean over measured-scale Y
            
        else;
            
            fprintf(1,'A_EstRtrAss: CndValM above UEP\n');
            EYM(iBS)=NaN;
            
        end;
        
    end; % bootstrap
    
    %Manage Inf
    EYM(isinf(EYM))=NaN;
    
    %We are interested in mean and median over bootstraps
    VlsEst.Ass(k,1)=mean(EYM(~isnan(EYM)));
    VlsEst.Ass(k,2)=median(EYM(~isnan(EYM)));
    
    %Record problematic bootstraps
    VlsEst.AssOK(k)=sum(~isnan(EYM));
    
end; % choice of return value


%% Estimate associated values 3-4 (ie q_61 and q_62), where the X return value is a function of bootstrap
% only do this for k=2
clear YL YM;
EYM=nan(nBS,1);

for iBS=1:nBS; %Loop over bootstrap resamples
    
    % estimate return value per bootstrap
    CndValL=VlsTru.Rtr.L;
    Prmh=Mdl.HT.Prmh(:,iBS);
    Rsdh=Mdl.HT.Rsd{iBS};
    %
    YL=A_YLSmp(CndValL,Prmh,Rsdh);
    %
    YM=pTrnScl('L2M',YL,[Mdl.GP.Prmh(:,2,iBS);0]);
    EYM(iBS)=mean(YM); %mean over measured-scale Y
    
end; % bootstrap

%Manage Inf
if sum(isinf(EYM))+sum(isnan(EYM))>0;
    fprintf('INF or NAN found!\n');
end;
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