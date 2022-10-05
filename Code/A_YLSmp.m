function YL=A_YLSmp(L,Prm,Rsd);

YL=Prm(1)*L+(L.^Prm(2))*(Prm(3) + Prm(4)*Rsd);

return;