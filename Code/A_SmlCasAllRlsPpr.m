function FrcBis=A_SmlCasAllRlsPpr(Dsg,Out,Cas);
% function FrcBis=A_SmlCasAllRls(Dsg,Out,Cas);
%
% Associated values: execute all realisations for a given case, and summarise results
%
% P. Jonathan, R. Towe 2022

%% Define case
Xix=Dsg.Xi(Cas(1));
Xiy=Dsg.Xi(Cas(2));
Dpn=Dsg.DpnNam{Cas(3)};
Kpp=Dsg.Dpn.(Dpn)(Cas(4));
SmpSiz=Dsg.SmpSiz(Cas(5));
HTThrU=Dsg.HTThrU(Cas(6));
AnnRat=Dsg.AnnRat; % annual rate of occurrence
RtrPrd=Dsg.RtrPrd; % return period of interest

%% Estimate true return values and associated values
Vls.Tru=A_ClcTru(Xix,Xiy,Dpn,Kpp,AnnRat,RtrPrd);

%% Loop over realisations
nMdl=1;
FB.Est=cell(nMdl,1); %first dimension MUST be the number of models
for iRls=1:Dsg.nRls;
    
    iM=0;
    
    %% Set counter for loop over model estimation schemes
    fprintf(1,'Cas %g%g%g%g%g%g Rls%g\n',Cas,iRls);
    
    %% Generate data on uniform margins from appropriate copula
    Dat=A_MakDat(SmpSiz,Dpn,Kpp,Xix,Xiy);
    
    %% Estimate HT model (HT1 uses Phil's fmincon)
    if 1;
        iM=iM+1; % Next model
        Mdl.Est{iM,1}=A_EstHT1(Dat,HTThrU,Dsg.nBS);
        Vls.Est{iM,1}=A_EstRtrAss(Mdl.Est{iM},Vls.Tru,AnnRat,RtrPrd);
    end;

    %% Estimate model using correct copula form with estimated margins
    if 1;
        iM=iM+1;
        Mdl.Est{iM,1}=A_EstCpl(Dat,Dpn,Dsg.nBS);
        Vls.Est{iM,1}=A_EstRtrAssCpl(Mdl.Est{iM},Dpn,AnnRat,RtrPrd);
    end;
    
    %% Estimate model using correct copula form with true margins
    if 1;
        iM=iM+1;
        Mdl.Est{iM,1}=A_EstCplTruMrg(Dat,Dpn,Dsg.nBS,Xix,Xiy);
        Vls.Est{iM,1}=A_EstRtrAssCpl(Mdl.Est{iM},Dpn,AnnRat,RtrPrd);
    end;
    
    %% Estimate HT model (HT1 uses Phil's fmincon) with true margins
    if 1;
        iM=iM+1; % Next model
        Mdl.Est{iM,1}=A_EstHT1TruMrg(Dat,HTThrU,Dsg.nBS,Xix,Xiy);
        Vls.Est{iM,1}=A_EstRtrAss(Mdl.Est{iM},Vls.Tru,AnnRat,RtrPrd);
    end;

    %% Estimate model using linear regression on physical scale with dependence threshold
    if 1;
        iM=iM+1;
        Mdl.Est{iM,1}=A_EstRgrThr(Dat,Dsg.nBS,HTThrU);
        Vls.Est{iM,1}=A_EstRtrAssRgr(Mdl.Est{iM},AnnRat,RtrPrd);
    end;

    %% Estimate bias
    FB=A_EstFrcBis(Vls,FB,iRls);
    
end; %iRls

%% Estimate bias and standard deviation statistics
FrcBis.Est=A_SmmFrcBis(FB);
FrcBis.Out=Out;
FrcBis.Cas=Cas;

%% Test case if no input arguments
if nargin==0;
    clf;
    subplot(1,2,1);
    imagesc(FrcBis.Est{1}.Ass.Raw(:,[7 19 8 20])'); colorbar;
    subplot(1,2,2);
    imagesc(FrcBis.Est{1}.AssOK.Raw(:,6)'); colorbar;
    drawnow;
end;

%% Save results
tStr=sprintf('FrcBis_Xix%g_Xiy%g_Dpn%s_Tht%g_SmpSiz%g_HTThrU%g.mat',Xix, Xiy, Dpn, Kpp, SmpSiz, HTThrU);
save(fullfile(Out,tStr),'FrcBis');

return;