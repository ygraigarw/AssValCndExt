function VlsTru=A_ClcTru(Xix,Xiy,Dpn,Kpp,AnnRat,RtrPrd);
%function VlsTru=A_ClcTru(Xix,Xiy,Dpn,Kpp,AnnRat,RtrPrd);
%
% Associated values: calculate true return and associated values
%
% P. Jonathan, R. Towe 2022

%% True return value for X on U, L and M scales
VlsTru.Rtr.U=1-(-1/AnnRat)*log(1-1/RtrPrd); %know this in closed form
VlsTru.Rtr.L=pTrnScl('U2L',VlsTru.Rtr.U);
VlsTru.Rtr.M=pTrnScl('U2M',VlsTru.Rtr.U,[Xix;1;0]);

%% True associated values on M scale using definitions 1 and 2
switch Dpn;
    case 'Gss'; % normal
        nRls=1e6;
        XNCnd=pTrnScl('U2N',VlsTru.Rtr.U);
        Rho=Kpp;
        PrmGP=[Xiy;1];
        Smm=A_SmmYMgXN(nRls,XNCnd,Rho,PrmGP);
        VlsTru.Ass.M(1)=Smm.Sml.Exp;
        VlsTru.Ass.M(2)=Smm.Sml.Mdn;
    case 'Lgs'; % logistic / Gumbel
        nRls=1e4;
        PrmGP=[Xiy;1];
        Alp=Kpp;
        XFCnd=pTrnScl('U2F',VlsTru.Rtr.U,PrmGP);
        Smm=A_SmmYMgXF(nRls,XFCnd,Alp,PrmGP);
        VlsTru.Ass.M(1)=Smm.Int.Exp;
        VlsTru.Ass.M(2)=Smm.Int.Mdn;
end;

%% Complete
return;