function FrcBis=A_SmmFrcBis(FB);
% function A_SmlCasAllRls(Dsg,Out,Cas);
%
% Associated values: summarise fractional bias
%
% P. Jonathan, R. Towe 2022

nMdl=size(FB.Est,1);
nRls=size(FB.Est{1}.Rtr,2);
FrcBis=cell(nMdl,1);
for iM=1:nMdl;
    
    %Bias in return value over all realisations
    FrcBis{iM}.Rtr.Men=mean(FB.Est{iM}.Rtr,2);
    
    %Bias in associated value over all realisations
    FrcBis{iM}.Ass.Men=mean(FB.Est{iM}.Ass,3);
    
    %Standard deviation in return value over all realisations
    FrcBis{iM}.Rtr.StnDvt=std(FB.Est{iM}.Rtr,[],2);
    
    %Standard deviation in associated value over all realisations
    FrcBis{iM}.Ass.StnDvt=std(FB.Est{iM}.Ass,[],3);
    
    %Return values per realisation
    FrcBis{iM}.Rtr.Raw=FB.Est{iM}.Rtr';
    
    %Associated values per realisation
    FrcBis{iM}.Ass.Raw=reshape(permute(FB.Est{iM}.Ass,[3,2,1]),nRls,20,1); 
    
    %Return values per realisation
    FrcBis{iM}.AssOK.Raw=FB.Est{iM}.AssOK'; 
    
end;

%% Complete
return;