clear; clc; clf;

%% GoPltBoxWhs
[jnk,tStr]=system('hostname'); CmpNam=tStr(1:end-1);
if strcmp(CmpNam,'LONBS-L-78932'); 
    Out='C:\Philip\Gwaith-Extremes\2022_MetoceanPaperWriting\MetOcean_Papers\AssociatedValues\SimulationStudy\Output';
    Smm='C:\Philip\Gwaith-Extremes\2022_MetoceanPaperWriting\MetOcean_Papers\AssociatedValues\SimulationStudy\Summary';
else;
    Out='Z:\project\MetOcean\AssociatedValues\MetOcean_Papers\AssociatedValues\SimulationStudy\Output_1000x100Ppr';
    Smm='Z:\project\MetOcean\AssociatedValues\MetOcean_Papers\AssociatedValues\SimulationStudy\Summary';
    addpath('C:\PhilipGit\pMtlUtl'); 
end;

if 1;
    load(fullfile(Out,'SmmFrcBis.mat'),'FrcBis','Dsg','Cas','Fal');
else;
    A=dir(Out);
    jA=0;
    for iA=1:size(A,1);
        if isempty(strfind(A(iA).name,'FrcBis_'))==0;
            jA=jA+1;
            load(fullfile(Out,A(iA).name));
            FrcBis0{jA,1}=FrcBis;
            Cas(jA,:)=FrcBis.Cas;
            clear FB;
        end;
    end;
    clear FrcBis;
    FrcBis=FrcBis0;
    Fal=zeros(size(FrcBis,1),1);
end;
nCas=size(FrcBis,1);

%% Case definition is: [iXix iXiy iDpn iTht iSmp iHT]%
I=ones(nCas,1);

%% Grand summary of performance of different estimators

% Superimpose multiple models on same figure
clf;
nMdl=size(FrcBis{1}.Est,1);
for iM=1:nMdl;
    
    % Return values
    Men.RV=nan(sum(I),5);
    for k=1:5;
        jC=0;
        for iC=1:nCas;
            if I(iC)==1 && Fal(iC)==0;
                nRls=size(FrcBis{iC}.Est{iM}.Rtr.Raw,1);
                jC=jC+1;
                Men.RV(jC:jC+nRls-1,k)=FrcBis{iC}.Est{iM}.Rtr.Raw(:,k);
                if FrcBis{iC}.Est{iM}.Rtr.Raw(:,k)>10;
                    fprintf(1,'Hello\n');
                end;
                jC=jC+nRls-1;
            end;
        end;
    end;
    
    Men.Ass=nan(sum(I),20);
    for k=1:20;
        jC=0;
        for iC=1:nCas;
            if I(iC)==1 && Fal(iC)==0;
                nRls=size(FrcBis{iC}.Est{iM}.Rtr.Raw,1);
                jC=jC+1;
                Men.Ass(jC:jC+nRls-1,k)=FrcBis{iC}.Est{iM}.Ass.Raw(:,k);
                jC=jC+nRls-1;
            end;
        end;
    end;
    
    Dat.Prm.Wdt=0.1;
    Dat.Prm.Thc=1;
    Dat.Prm.YLmt=[-2 2];
    Dat.Prm.Clr=pClr(iM);
    Dat.Prm.Off=2*(iM-(nMdl+1)/2)/10;
    Dat.Rtr.Err=Men.RV;
    Dat.Ass.Err=Men.Ass;
    A_BoxWhs(Dat);
    
end;

pDatStm('GrandSummary');
pGI(fullfile(Smm,'Smm1GrandSummary'),2);