%% USER NOTE on PRO+
%Code to calculate protrusion and critical shear stresses for a distribution of grains on the bed.
%Code must be run after G3Point (Steer et al., 2022) and requires variables from
%G3Point as well as the G3Point modified point cloud. 
%Use this code at your own risk and modifications may be needed for your
%specific application. If you encounter errors in the code please tell
%Elowyn Yager (eyager@uidaho.edu). Code written by Elowyn Yager and last
%modified on 11/15/2023.
%% Read input variables and inputs from G3Point/other software, calculate protrusion 
clearvars
[proinputs]=defineinputs; %read inputs from .csv file
inputname=string(proinputs.inputfile);savepro=append('protrusion',inputname);    

if proinputs.whichcalcs==1%Determine whether to step through all calcs or just the critical shear
    %stress calcs (only possible if you have previously completed protrusion
    %calcs for same data)
    
%Load input varibles from G3Point.other software    
    load(inputname); 
% Calculate point cloud and grain size metrics
    if proinputs.whichinput==1 %if input is from G3Point
        b=double(granulo.diameter(2,:))'; %entire GSD that correponds to each protrusion value
    end
    b(b<0.002)=NaN; %eliminate any sand or finer particles from the distribution
    b84=prctile(b,84); %84th percentile of the GSD
    b50=prctile(b,50); %median of the GSD 
    bedstd=std(ptCloud.Location(:,3)); %standard deviation of bed elevation from point cloud
    
% Calculate driving (pD, 10th percentile) and resisting (pR, median) protrusion for each grain on the bed
    if proinputs.whichinput==1 %if input is from G3Point
    [pR,pD]=protrusionG3calc(proinputs,bedstd,b84,ptCloud,labels,nlabels,Ellipsoidm); 
    else
    [pR,pD]=protrusionothercalc(proinputs,bedstd,b84,ptCloud,grain); 
    end
    
%Save the protrusion results such that does not need to be calculated again
    %if want to vary critical shear stress variables. 
    save(savepro,"pR", "pD", "b", "b84", "b50", "bedstd")
    
elseif proinputs.whichcalcs==0
    load(savepro); 
end
%% Calculate driving and resisting forces and critical shear stress for each grain on the bed
%details on these function outputs are provided in taucrcalc.m. Some
%variables output by taucrcalc.m are not used here and are therefore not
%specified.
[FrWt_dist,taucr_dist,taustcr_dist,~,~,b_dist,pD_dist,pR_dist]=taucrcalc(proinputs,b,b50,b84,bedstd,pD,pR);
%% Calculate critical shear stress and protrusion for each grain size bin
[Binned,bbinmean,bbins]=bintau(b,pD,pR,b_dist,pD_dist,pR_dist,taucr_dist,taustcr_dist,proinputs);
%% Ouput calculated variables
if proinputs.saveoutputs==1  
    outputname=string(proinputs.outputfile);
    save(outputname, "Binned", "bbins", "bbinmean", "b_dist", "FrWt_dist", "taucr_dist", "taustcr_dist", "b", "pD", "pR","proinputs")
end

