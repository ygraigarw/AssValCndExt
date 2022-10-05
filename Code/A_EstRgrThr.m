function Mdl=A_EstRgrThr(Dat,nBS,HTThrU);

SmpSiz=size(Dat.M,1);

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
    
    %% Marginal model for X only
    tGPPrmX=gpfit(YBS(:,1))';
    %Add to structure
    Mdl.GP.Prmh(:,1,iBS)=tGPPrmX;
    
    %% Use only data above threshold for linear regression
    tThr=quantile(YBS(:,1),HTThrU);
    tYBS=YBS(YBS(:,1)>tThr,:);
    tSmpSiz=size(tYBS,1);

    %% Linear regression
    Prmh=[ones(tSmpSiz,1) tYBS(:,1)]\tYBS(:,2);  
    %Add to structure (save the intercept and slope terms)
    Mdl.Cpl.Prmh(:,:,iBS)=Prmh;
    
end;

return;