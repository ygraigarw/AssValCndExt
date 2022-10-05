function Smm=A_SmmYMgXF(nRls,XFCnd,Alp,PrmGP);
%function Smm=A_SmmYMgXF(nRls,XFCnd,Alp,PrmGP);
%
% Conditional mean and median on Measurement scale for bivariate Logistic dependence with parameter Alp.
% Marginal distribution for bivariate Logistic is Frechet.
% Measurement scale is GP with parameter PrmGP.
% Conditioning value on Frechet scale is XNCnd.
%
% The approach used is motivated by Joe p172 S4.8.1.
%
% P. Jonathan, R. Towe 2022

%% Test case
if nargin==0;
    %nRls=1e3;%1e8;
    %XFCnd=20;
    XFCnd=-1/log(1-1e-7);
    Alp=0.9;
    PrmGP=[-0.2;1];
end;

%% Calculate conditional expectation and median in closed form (found to perform best overall)
if 1;
    [Smm.Int.Exp,Smm.Int.Mdn]=A_EYMgXF(XFCnd,Alp,PrmGP);
end;

%% Simulate data from conditional bivariate logistic (not used for paper)
if 0; 
    
    %Expression h=0 to solve for each choice of random u at conditioning value x
    %Solution gives point from conditional distribution for bivariate logistic (or bivariate Gumbel)
    h=@(v,Alp,x,u) v+(1/Alp-1)*log(v)+(1/Alp-1)*log(x)-1/x+log(u);
    dh=@(v,Alp) 1+(1/Alp-1)./v;
    
    %Newton-Raphson solve to generate each observation
    V=nan(nRls,1);
    for i=1:nRls;
        
        U=rand;
        
        %Find valid starting solutions
        vGrd=linspace(0,10,10000)'; %Grid
        DltGrd=h(vGrd,Alp,XFCnd,U)./dh(vGrd,Alp); %Newton-Raphson step on the grid of starting options
        GrdStr=vGrd(imag(DltGrd)==0); %Valid steps
        hStr=h(GrdStr,Alp,XFCnd,U);
        if nargin==0;
            subplot(1,3,1);plot(vGrd,h(vGrd,Alp,XFCnd,U),'k');
            subplot(1,3,2);plot(vGrd,dh(vGrd,Alp),'k');
            subplot(1,3,3);plot(GrdStr,hStr);
        end;
        
        %Newton-Raphson
        for j=1:nItr;
            
            if j==1; %Sensible valid starting solution
                t=find(abs(hStr)==min(abs(hStr)));t=t(1);
                tV=GrdStr(t);
            else; %Newton-Raphson iteration
                Dlt=h(tV,Alp,XFCnd,U)/dh(tV,Alp);
                if imag(Dlt)~=0;
                    fprintf(1,'Error - imaginary solution\n');
                end;
                if abs(Dlt)>1e-12; %Solution found
                    if tV>Dlt;
                        tV=tV-Dlt;
                    else;
                        tV=tV/10;
                    end;
                else;
                    break;
                end;
            end;
            
        end;
        if j==nItr; %Invalid solution since max iterations
            fprintf(1,'Warning: max iters\n');
            V(i)=nan;
        else;
            V(i)=tV;
        end;
        
        if 0; %Counter to StdOut
            if rem(i,1000)==0;
                fprintf(1,'+');
            end;
            if rem(i,10000)==0;
                fprintf(1,'\n');
            end;
        end;
        
    end;
    
    %Solve for Fr (given V and X)
    t=V.^(1/Alp)-XFCnd.^(-1/Alp);
    if sum(t<0)>0; %Check for bad values (t is theoretically >0, just numerical issues cause this)
        fprintf(1,'Error - negative value\n');
        t(t<0)=nan;
    end;
    XF=(t).^(-Alp);
    
    %Output simulated mean and median on measured scale
    XM=pTrnScl('F2M',XF,PrmGP);
    Smm.Sml.Exp=mean(XM);
    Smm.Sml.Mdn=median(XM);
    
    %Output simulated mean and median on Frechet scale
    fprintf(1,'Conditioning value on measurement scale=%g\n',pTrnScl('F2M',XFCnd,PrmGP));
    fprintf(1,'On measured scale: Sml(Exp=%g,Mdn=%g)\n',mean(XM),median(XM));
    fprintf(1,'On measured scale: Int(Exp=%g,Mdn=%g)\n',Smm.Int.Exp,Smm.Int.Mdn);
    
end;

%% Complete
return;
