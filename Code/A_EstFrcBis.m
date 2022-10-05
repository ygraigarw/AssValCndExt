function FB=A_EstFrcBis(Vls,FB,iRls);
% function A_SmlCasAllRls(Dsg,Out,Cas);
%
% Associated values: estimate fractional bias
%
% P. Jonathan, R. Towe 2022

for iM=1:size(Vls.Est,1);
    for k=1:5;
        FB.Est{iM,1}.Rtr(k,iRls)=Vls.Est{iM}.Rtr(k)./Vls.Tru.Rtr.M-1;
        %You are interested in comparing with E_Y (estimated using E_Z(E_Y)
        %or med_Z(E_Y). You are not interested in med_Y!!!!!!
        tAssTru=ones(4,1)*Vls.Tru.Ass.M(1);
        FB.Est{iM,1}.Ass(k,:,iRls)=Vls.Est{iM}.Ass(k,:)./tAssTru'-1;

    end;
    % Save record of failed runs
    FB.Est{iM,1}.AssOK(:,iRls)=Vls.Est{iM}.AssOK';
end;

%% Complete
return;