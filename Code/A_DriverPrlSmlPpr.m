% ------------------------------------------------------
%% A_DriverPrlSml
%
% Associated values: simulation driver
%
% P. Jonathan, R. Towe 2022
%
% Driver currently set up for quick test run.
% To run full simulation, need to switch parallel on (IsPrl=1) and use full
% design (see lines 52-71)
%
% Warning: full simulation takes approximately 1-2 weeks running 40 cores
% in parallel.

%% Simulation of data
% MARGINS
% Grid of marginal xi for X      -0.4:0.2:0.2
% Grid of marginal xi for Y      -0.4:0.2:0.2
% Marginal GP threshold = 0 for X and Y
% Marginal GP scale = 1 for X and Y
% DEPENDENCE TYPE, D
% Gaussian             "rho"       0.1 0.5 0.9 (kappa in paper)
% Logistic             "alpha"     0.9 0.5 0.1 (1-kappa in paper)
% SAMPLE SIZE
% n=200, 1000, 10k
% REALISATIONS
% nRls=1000
% COMPUTATIONAL COST
% 4 xiX * 4 xiY * 2 D * 3 kappa * 3 n * 1000 nRls

%% Model estimation
% BOOTSTRAPS
% nBS=100
% THRESHOLD NEPs
% (Marginal=0)
% Conditional   tau=0.5, 0.8, 0.9

%% Total computational cost
% 4 xiX * 4 xiY * 2 D * 3 kappa * 3 n * 3 tau * 1000 nRls * 100 nBS
% ------------------------------------------------------

%% Set-up
% Clear everything
clf; clear; clc;

%% Folder for output
Out='..\SimulationOutput'; % Folder for output data % SET THIS AS DESIRED

%% Control parallel execution
% Use parallel processing for full design, but serial run OK for checking
nWrk=40; % Number of workers
IsPrl=0; % Specify whether to run in parallel (=1) or serial (=0)

%% Set up experimental design structure
if 1; %Quick run for checking (THIS WILL TAKE ABOUT A MINUTE WITH A SINGLE CORE)
    Dsg.Xi=[-0.4 0.2]';%(-0.4:0.2:0.2)';
    Dsg.DpnNam={'Gss';'Lgs'};
    Dsg.Dpn.Gss=[0.1 0.9]';%[0.1 0.5 0.9]';
    Dsg.Dpn.Lgs=[0.9 0.1]';%[0.9 0.5 0.1]';
    Dsg.SmpSiz=1000; %[200 1000 10000]';
    Dsg.HTThrU=0.9;%[0.5 0.8 0.9]';
    Dsg.nRls=5;
    Dsg.nBS=5;
end;
if 0; %Full run (THIS WILL TAKE OF THE ORDER OF 2 WEEKS WITH 40 CORES)
    Dsg.Xi=(-0.4:0.2:0.2)';
    Dsg.DpnNam={'Gss';'Lgs'};
    Dsg.Dpn.Gss=[0.1 0.5 0.9]';
    Dsg.Dpn.Lgs=[0.9 0.5 0.1]';
    Dsg.SmpSiz=[200 1000 10000]';
    Dsg.HTThrU=[0.5 0.8 0.9]';
    Dsg.nRls=1000;
    Dsg.nBS=100;
end;

% Set up return anad associate values structure
Dsg.AnnRat=100;
Dsg.RtrPrd=100;

% save Dsg
save(fullfile(Out,'Dsg.mat'),'Dsg');

%% Set up case indicator over design
iCas=0;
nCas=size(Dsg.Xi,1)*size(Dsg.Xi,1)*size(Dsg.DpnNam,1)*size(Dsg.Dpn.Gss,1)*size(Dsg.SmpSiz,1)*size(Dsg.HTThrU,1);
Cas=nan(nCas,6);
for iXix=1:size(Dsg.Xi,1);
    for iXiy=1:size(Dsg.Xi,1);
        for iDpn=1:size(Dsg.DpnNam,1);
            for iTht=1:size(Dsg.Dpn.Gss,1); % Assumes that there are common numbers of theta for Gss and Lgs
                for iSmp=1:size(Dsg.SmpSiz,1);
                    for iHT=1:size(Dsg.HTThrU,1);
                        iCas=iCas+1;
                        Cas(iCas,:)=[iXix iXiy iDpn iTht iSmp iHT];
                    end;
                end;
            end;
        end;
    end;
end;

%% Run simulation
Fal=nan(nCas,1);
Err=cell(nCas,1);
if IsPrl==0; %Series (DO NOT USE PARALLEL PROCESSING)
    
    % Serial simulation
    FrcBis=cell(nCas,1);
    for iCas=1:nCas;
        FrcBis{iCas}=A_SmlCasAllRlsPpr(Dsg,Out,Cas(iCas,:));
    end;
    Fal=zeros(nCas,1);
    
else; %Parallel
    
    if 1; % New parallel pool
        delete(gcp('nocreate')); %delete any existing pool
        h=parcluster; %create a new cluster
        h.NumWorkers=nWrk; %specify number of workers on laptop
        parpool(h.NumWorkers); %create new parallel pool
    end;
    
    % Parallel simulation
    FrcBis=cell(nCas,1);
    parfor iCas=1:nCas; %Parfor loop, run and save each case independently
        
        % Check if results already present; if so use them
        tCas=Cas(iCas,:);
        Xix=Dsg.Xi(tCas(1));
        Xiy=Dsg.Xi(tCas(2));
        Dpn=Dsg.DpnNam{tCas(3)};
        Tht=Dsg.Dpn.(Dpn)(tCas(4));
        SmpSiz=Dsg.SmpSiz(tCas(5));
        HTThrU=Dsg.HTThrU(tCas(6));

        tStr=sprintf('FrcBis_Xix%g_Xiy%g_Dpn%s_Tht%g_SmpSiz%g_HTThrU%g.mat',Xix, Xiy, Dpn, Tht, SmpSiz, HTThrU);
        
        % Check if output file exists; if so, don't overwrite
        if exist(fullfile(Out,tStr),'file')==0; 
            try;
                FrcBis{iCas}=A_SmlCasAllRlsPpr(Dsg,Out,Cas(iCas,:));
                Fal(iCas)=0;
                fprintf(1,'Successful case %g%g%g%g%g%g\n',Cas(iCas,:));
            catch ME;
                FrcBis{iCas}.Err=ME;
                Fal(iCas)=1;
                fprintf(1,'*****Failure for case %g%g%g%g%g%g*****\n',Cas(iCas,:));
            end;
        else; %Load existing results
            tFrcBis=load(fullfile(Out,tStr),'FrcBis');
            FrcBis{iCas}=tFrcBis.FrcBis;
            Fal(iCas)=-1;
            fprintf(1,'*****Prior result for case %g%g%g%g%g%g*****\n',Cas(iCas,:));
        end;
        
    end;
    
end;

%% Save full results in one MATLAB file
save(fullfile(Out,'SmmFrcBis.mat'),'FrcBis','Dsg','Cas','Fal');