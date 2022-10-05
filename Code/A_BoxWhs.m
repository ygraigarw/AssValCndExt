function A_BoxWhs(Dat);

if nargin==0;
    Dat.Prm.Wdt=0.2;
    Dat.Prm.Thc=1;
    Dat.Prm.YLmt=[-0.5 0.5];
    Dat.Prm.Clr=[1 0 0];
    Dat.Prm.Off=0.1;
    Dat.Rtr.Err=randn(100,5);
    Dat.Ass.Err=randn(100,20); Dat.Ass.Err(:,[3 4 11 12 15 16])=nan;
    Dat.Txt='Test';
end;

pLtx;

%% Box-whisker return value
subplot(2,1,1);hold on;
set(gca,'xlim',[1 5]);
xlmt=get(gca,'xlim');%plot(xlmt,[0 0],'-','color',0.1*ones(3,1),'linewidth',1,'linestyle','--');
Lct=0;
for k=1:5;
    Lct=Lct+1;
    if sum(isnan(Dat.Rtr.Err(:,Lct))==0)>0;
        pBoxWhs(Lct+Dat.Prm.Off, Dat.Prm.Wdt, Dat.Rtr.Err(:,k), Dat.Prm.Clr, Dat.Prm.Thc);
    end;
end;
tTck={'$q_1$';'$q_2$';'$q_3$';'$q_4$';'$q_5$';};
set(gca,'xtick',(1:5),'xticklabel',tTck);
pAxsLmt; pDflBig;
ylabel 'FB return value';
ylmt=get(gca,'ylim');
set(gca,'ylim',[max(ylmt(1),Dat.Prm.YLmt(1)) min(ylmt(2),Dat.Prm.YLmt(2))]);

%% Box-whisker associated value
subplot(2,1,2);hold on;
set(gca,'xlim',[1 20]);
%xlmt=get(gca,'xlim');%plot(xlmt,[0 0],'-','color',0.1*ones(3,1),'linewidth',1,'linestyle','--');
Lct=0;
Lct2=0;
for k=1:5;
    for j=1:4;
        Lct=Lct+1;
        if sum(isnan(Dat.Ass.Err(:,Lct))==0)>0;
            Lct2=Lct2+1;
            pBoxWhs(Lct2+Dat.Prm.Off, Dat.Prm.Wdt, Dat.Ass.Err(:,Lct), Dat.Prm.Clr, Dat.Prm.Thc);
        end;
    end;
end;
tTck={'$q_{11}^*$';'$q_{12}^*$';...
    '$q_{21}^*$';'$q_{22}^*$';'$q_{23}^*$';'$q_{24}^*$';...
    '$q_{31}^*$';'$q_{32}^*$';...
    '$q_{41}^*$';'$q_{42}^*$';...
    '$q_{51}^*$';'$q_{52}^*$';'$q_{53}^*$';'$q_{54}^*$';};
set(gca,'xtick',(1:14),'xticklabel',tTck);
pAxsLmt; pDflBig;
ylabel('FB associated value');
ylmt=get(gca,'ylim');
set(gca,'ylim',[max(ylmt(1),Dat.Prm.YLmt(1)) min(ylmt(2),Dat.Prm.YLmt(2))]);

return;