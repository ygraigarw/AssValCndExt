function CndDns=A_fYFgXF(YF,XFCnd,Alp);
% function CndDns=A_fYFgXF(YF,XFCnd,Alp);
%
% Conditional density of the bivariate logistic (also called bivariate Gumbel)
% Variables on Frechet scale by definition
%
% P. Jonathan, R. Towe 2022

%% Test case
if nargin==0; % Test
    XFCnd=1e20;
    DltF=XFCnd*1e-3;
    YFMxm=XFCnd*1e2;
    YF=(DltF:DltF:YFMxm)';
    Alp=0.08;
end;

%% Basic calculation on Frechet scale
% Better to work with log conditional density, then
% exponentiate at the end
% Idea is:
% G=(XFCnd.^(-1/Alp)+YF.^(-1/Alp));
% V=G.^Alp; % Exponent measure
% FXY=exp(-V); % Joint distribution
% dGdx=(-1/Alp)*XFCnd.^(-1/Alp-1);
% dGdy=(-1/Alp)*YF.^(-1/Alp-1);
% dVdx=Alp*G.^(Alp-1).*dGdx;
% dVdy=Alp*G.^(Alp-1).*dGdy;
% d2Vdxdy=Alp*(Alp-1)*G.^(Alp-2).*dGdx.*dGdy;
% fXY=FXY.*(dVdx.*dVdy-d2Vdxdy);
% CndDns=fXY.*(XFCnd.^2).*exp(1./XFCnd);

% G=(XFCnd.^(-1/Alp)+YF.^(-1/Alp));
% LgrCndDns = ...
%     - G.^Alp ...
%     + (Alp-2).*log(G) ...
%     - (1/Alp - 1).*log(XFCnd) ...
%     - (1/Alp + 1).*log(YF) ...
%     + log(G.^Alp+(1-Alp)/Alp) ...
%     + 1/XFCnd;

G=(XFCnd.^(-1/Alp)+YF.^(-1/Alp));
GPA=G.^Alp;
LgrCndDns = ...
    - GPA ...
    + (Alp-2).*log(G) ...
    - (1/Alp - 1).*log(XFCnd) ...
    - (1/Alp + 1).*log(YF) ...
    + log(GPA+(1-Alp)/Alp) ...
    + 1/XFCnd;

CndDns=exp(LgrCndDns);

%% Test case
if nargin==0;
    clf; hold on;
    %plot(YF,CndDns,'k-');
    plot(YF,exp(LgrCndDns),'r--');
end;

%% Complete
return;
