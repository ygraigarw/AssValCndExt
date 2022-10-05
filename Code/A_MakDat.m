function Dat=A_MakDat(n,Dpn,Kpp,Xix,Xiy);
%function Dat=A_MakDat(n,Dpn,Kpp,Xix,Xiy);
%
% Associated values: simulate data
%
% P. Jonathan, R. Towe 2022

%% Generate data on uniform margins from appropriate copula
switch Dpn;
    case 'Gss';
        rho=Kpp;
        Z1=randn(n,1);
        Z2=randn(n,1).*sqrt(1-rho.^2)+rho.*Z1;
        Z=[Z1 Z2];
        Dat.U=normcdf(Z); %uniform;
    case 'Lgs';
        alp=Kpp;
        Dat.F=A_RndSmlLgs(alp,2,n)'; %function assumes Frechet margins
        Dat.U=pTrnScl('F2U',Dat.F);
end;
Dat.L=pTrnScl('U2L',Dat.U);

%% PIT to generate marginal values
Dat.M(:,1)=pTrnScl('U2M',Dat.U(:,1),[Xix;1;0]);
Dat.M(:,2)=pTrnScl('U2M',Dat.U(:,2),[Xiy;1;0]);

%% Complete
return;