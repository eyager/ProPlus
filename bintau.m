function [Binned, bbinmean,bbins]=bintau(b,pD,pR,b_dist,pD_dist,pR_dist,taucr_dist,taustcr_dist,proinputs)
%% USER NOTE
%Code to take estimated critical shear stresses from taucrcalc and estimated
%protrusion values from protrusioncalc and dividethem into 
%different grain size bins. Critical Shields stresses are also used to estimate
%a hiding function
%Use this code at your own risk and modifications may be needed for your
%specific application. If you encounter errors in the code please tell
%Elowyn Yager (eyager@uidaho.edu).  Code written by Elowyn Yager and last
%modified on 11/15/2023.
%% SET UP GRAIN SIZE BINS AND STRUCTURE FOR VARIABLES FOR EACH GRAIN SIZE BIN
disp('--- BINNING CRITICAL SHEAR STRESSES AND PROTRUSION BY GRAIN SIZE');
phi=-1:proinputs.phistep:-11; %grain size bins determined for 
bbins=(2.^(-phi))./1000; %grain size bins in units of meters. Bins range 
    %in size from 0.002 m (lowest size recommended to use in Pro+) to 2.048
    %m. If larger boulders are desired in the distribution, the phi range
    %can be changed. 
bbinmean=geomean([bbins(1:end-1); bbins(2:end)]); %geometric mean of each grain size bin    
    
%eliminate critical shear stresses and cooresponding protrusion values
%for which estimated critical shear stresses could not be calculated (i.e.
%resulted in a NaN for critical shear stress). This may produce a different
%protrusion distribution than that originally calculated and therefore two
%protrusion distributions for pDi and pRi (for each ith grain size) are provided: the original
%calculated distributions (see below) and the distributions that correspond to those
%critical shear stress values that were numbers (not NaN). pDi and pRi
%distributions associated with taucr are labeled pDi_for_taucr and
%pRi_for_taucr, respectively. 
indnan=isnan(taucr_dist);taucr_dist(indnan==1)=[];taustcr_dist(indnan==1)=[];
b_dist(indnan==1)=[];pD_dist(indnan==1)=[];pR_dist(indnan==1)=[];

%find indices of the expanded grain size distribution (from the critical shear stress
%calculations) that have grain sizes within the range of values allowed in 
%each grain size bin. These indices correspond to those needed for 
%taucr_dist and taustcr_dist in each grain size bin.  
indtaucr=discretize(b_dist,bbins);Nbinstaucr=histcounts(b_dist,bbins);

%find indices of the original grain size distribution that have grain sizes
%within the range of values allowed in each grain size bin. These indices
%correspond to those needed for original pD and pR distributions for each grain size bin.  
indp=discretize(b,bbins);Nbinsp=histcounts(b,bbins); 


%setup structure to house critical shear stress (taucri), critical Shields
%stress (taustcri), driving protrusion (pDi), and resisting protrusion
%(pRi), driving protrusion associated with each taucr value (pDi_for_taucri),
% and resisting protrusion associated with each taucr value (pRi_for_taucri), 
%for each grain size bin
svalue1={NaN.*ones(length(bbins)-1,max(Nbinsp))};
svalue2={NaN.*ones(length(bbins)-1,max(Nbinstaucr))};
Binned=struct('taucri',svalue2,'taustcri',svalue2,'pDi_for_taucri',svalue2,'pRi_for_taucri',svalue2,'pDi',svalue1,'pRi',svalue1);
clear svalue1 svalue2
%% BIN CRITICAL SHEAR STRESSES AND PROTRUSIONS BY GRAIN SIZE BIN
%binned critical shear stresses
for j=1:nanmax(indtaucr)
    Binned.taucri(j,1:numel(b_dist(indtaucr==j)))=taucr_dist(indtaucr==j);
    Binned.taustcri(j,1:numel(b_dist(indtaucr==j)))=taustcr_dist(indtaucr==j);
    Binned.pDi_for_taucri(j,1:numel(b_dist(indtaucr==j)))=pD_dist(indtaucr==j);
    Binned.pRi_for_taucri(j,1:numel(b_dist(indtaucr==j)))=pR_dist(indtaucr==j);
end
%binned protrusions
for j=1:nanmax(indp)
    Binned.pDi(j,1:numel(b(indp==j)))=pD(indp==j);
    Binned.pRi(j,1:numel(b(indp==j)))=pR(indp==j);
end
end