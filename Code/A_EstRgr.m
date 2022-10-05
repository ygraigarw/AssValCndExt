function Mdl=A_EstRgr(Dat,nBS);

SmpSiz=size(Dat.M,1);

for iBS=1:nBS; % loop over bootstraps
    
    % bootstrap indicator
    if iBS==1;
        BSInd=(1:SmpSiz)';
    else;
        BSInd=randsample((1:SmpSiz)',SmpSiz,1);
    end;
    % add to structure
    Mdl.BSInd(:,iBS)=BSInd;
    
    % bootstrap resample
    YBS=Dat.M(BSInd,:);
    
    % marginal model for X only
    tGPPrmX=gpfit(YBS(:,1))';
    % add to structure
    Mdl.GP.Prmh(:,1,iBS)=tGPPrmX;
    
    % linear regression
    Prmh=[ones(SmpSiz,1) YBS(:,1)]\YBS(:,2);
        
    % add to structure
    % save the intercept and slope terms
    Mdl.Cpl.Prmh(:,:,iBS)=Prmh;
    
end;

return;