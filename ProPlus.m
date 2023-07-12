%% USER NOTE on PRO+
%Code to calculate protrusion and critical shear stresses for a distribution of grains on the bed.
%Code must be run after G3Point (Steer et al., 2022) and requires variables from
%G3Point as well as the G3Point modified point cloud. 
%Use this code at your own risk and modifications may be needed for your
%specific application. If you encounter errors in the code please tell
%Elowyn Yager (eyager@uidaho.edu). Code written by Elowyn Yager and last
%modified on 4/5/2023.
%% Read input variables and inputs from G3Point/other software
clearvars
[proinputs]=defineinputs;
if proinputs.whichinput==1 %if input is from G3Point
    load('G3Pointinput.mat');
elseif proinputs.whichinput==2 %if input is from other software
    load('otherinput.mat'); 
else
    disp('error, need to enter defineinput choice between 1-2')
end
%%  Calculate point cloud and grain size metrics
if proinputs.whichinput==1 %if input is from G3Point
    b=double(granulo.diameter(2,:))'; %entire GSD that correponds to each protrusion value
end
b(b<0.002)=NaN; %eliminate any sand or finer particles from the distribution
b84=prctile(b,84); %84th percentile of the GSD
b50=prctile(b,50); %median of the GSD 
bedstd=std(ptCloud.Location(:,3)); %standard deviation of bed elevation from point cloud
%% Calculate protrusion for each grain on the bed
%calculate driving (pD, 10th percentile) and resisting (pR, median) protrusion 
if proinputs.whichinput==1 %if input is from G3Point
   [pR,pD]=protrusionG3calc(proinputs,bedstd,b84,ptCloud,labels,nlabels,Ellipsoidm); 
else
   [pR,pD]=protrusionothercalc(proinputs,bedstd,b84,ptCloud,grain); 
end
%% Calculate driving and resisting forces and critical shear stress for each grain on the bed
%details on these function outputs are provided in taucrcalc.m. Some
%variables output by taucrcalc.m are not used here and are therefore not
%specified.
[FrWt_dist,taucr_dist,taustcr_dist,~,~,b_dist,pD_dist,pR_dist]=taucrcalc(proinputs,b,b50,b84,bedstd,pD,pR);
%% Calculate critical shear stress and protrusion for each grain size bin
[Binned,bbins]=bintau(b,pD,pR,b_dist,pD_dist,pR_dist,taucr_dist,taustcr_dist,proinputs);
%% Ouput calculated variables
if proinputs.saveoutputs==1
    save ProPlusOutputs.mat Binned bbins FrWt_dist taucr_dist taustcr_dist b pD pR
end
