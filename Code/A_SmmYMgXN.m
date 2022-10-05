function Smm=A_SmmYMgXN(nRls,XNCnd,Rho,PrmGP);
% function SmmYMgXN(XNCnd,Rho,PrmGP);
%
% Conditional mean and median on Measurement scale for bivariate Normal dependence with correlation Rho.
% Measurement scale is GP with parameter PrmGP.
% Conditioning value on Normal scale is XNCnd.
%
% P. Jonathan, R. Towe 2022

if nargin==0; % Test
    nRls=1e3;
    Rho=0.999;
    XNCnd=500;
    PrmGP=[-0.2;1];
end;

%% Simulate data from conditional Gaussian (found to be most stable in extreme cases)
YN=randn(nRls,1).*sqrt(1-Rho.^2)+Rho.*XNCnd;
M=pTrnScl('N2M',YN,PrmGP); %Transform to measured scale
Smm.Sml.Exp=mean(M);
Smm.Sml.Mdn=median(M);

%% Calculate mean and median by numerical integration (not used for paper)
if 0;
    % Grid
    XNMxm=100;
    DltN=0.01;
    try;
        XN=(-XNMxm:DltN:XNMxm)';
        Smm.Int.Exp=sum(A_fYMgXN(XN,XNCnd,Rho,PrmGP)*DltN);
        % Median (this is easy since the conditional distribution is Gaussian and so symmetric, so median=mean)
        YN=Rho*XNCnd;
        Smm.Int.Mdn=pTrnScl('N2M',YN,PrmGP);
    catch;
        Smm.Int.Exp=NaN;
        Smm.Int.Mdn=NaN;
    end;

    % Use other functions to check result
    Smm.Int1.Exp=A_EYMgXN1(XNCnd,Rho,PrmGP);
    Smm.Int2.Exp=A_EYMgXN2(XNCnd,Rho,PrmGP);
    
    fprintf(1,'Conditioning value on measurement scale=%g\n',pTrnScl('N2M',XNCnd,PrmGP));
    fprintf(1,'Sml.Exp=%g Sml.Mdn=%g\n',Smm.Sml.Exp,Smm.Sml.Mdn);
    fprintf(1,'Int.Exp=%g Int.Mdn=%g\n',Smm.Int.Exp,Smm.Int.Mdn);
    fprintf(1,'Int1.Exp=%g\nInt2.Exp=%g\n',Smm.Int1.Exp,Smm.Int2.Exp);
end;

%% Complete
return;



