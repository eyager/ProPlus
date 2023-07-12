function [proinputs]=defineinputs
%Reads the proinputs.csv file that located in the same directory as
% point cloud. Default values are below and are used if this file is not
% located. Places all input variables into structure called proinputs.

%% If the .csv file exists, then read the input variables for Pro+ from the input file
if isfile('proinputs.csv')
    disp('--- READING PRO+ INPUT VARIABLES FROM INPUT FILES');
    [proinputs.saveoutputs,proinputs.usepivot,proinputs.whichinput,proinputs.whichks, proinputs.whichradius,...
        proinputs.setks,proinputs.setradius,proinputs.perimpoints,proinputs.meanphif,...
        proinputs.stdphif, proinputs.numphif,proinputs.phistep, proinputs.rhos,...
        proinputs.rho,proinputs.porosity,proinputs.Cd,proinputs.Cl,proinputs.Cv,...
        proinputs.K,proinputs.g]=...
        readvars('proinputs.csv');
        %readvars([param.ptCloudpathname 'proinputs.csv']);
    %INPUT VARIABLES WITH CHOICE OF YES (1) OR NO (0)
    %saveoutputs: Specifies whether to save the protrusion and critical 
        %shear stress outputs (=1) or not (=0
    %usepivot: Specifies whether to include the pivot angle (=1) or not (=0)
        %include pivot angle in critical shear stress calculations
    %INPUT VARIABLES WITH CHOICE OF TWO TO THREE DIFFERENT OPTIONS (1,2,3) 
    %whichinput: Specifies the method of inputing the necessary grain
        %variables for protrusion calculations. G3Point (=1) supplies grain
        %diameters, the detrended overall point cloud, and the detrended point
        %cloud for each grain, which must be saved to G3Pointinput.mat
        %Other input (=2) is any other user defined method
        %of obtaining the detrended point cloud, and the diameters and
        %perimeter coordinates of each grain on the bed. If 2 is chosen,
        %details on the formatting of the inputs in the otherinputs.mat
        %file are provided in the user manual. 
    %whichks: Specifies what ks value to use in critical shear stress
        %calculations: the standard deviation of bed elevations from the 
        %point cloud (=1) or the calculated 84th percentile of the grain 
        %size distribution b84 (=2)
    %whichradius: Specifies what search radius  to use in protrusion
        %calculations: the standard deviation of bed elevations from the 
        %point cloud (=1), the calculated 84th percentile of the grain 
        %size distribution b84 (=2), or a user specfied value (=3)
    %OPTIONAL INPUT VARIABLES WITH USER SET VALUES 
    %setradius: Specifies the user defined value for the search radius to
        %use when this option (3) is chosen for whichradius. Enter a number in
        %units of meters. There must be a number here, this number will not
        %be used if 1 or 2 is specified in whichradius.
    %setks: Specifies the user defined value for ks  to
        %use when this option (3) is chosen for whichks. Enter a number in
        %units of meters. There must be a number here, this number will not
        %be used if 1 or 2 is specified in whichks.
    %REQUIRED INPUT VARIABLES WITH USER SET VALUES     
    %perimpoints: Specifies the point increment on the grain perimeter to use
        %in the protrusion surrounding bed area search.  Smaller point
        %increments result in a search shape that more closely mimic the
        %grain perimeter but larger increments result in faster
        %calculations.  Recommend using the smallest increment that is
        %feasible. Increments greater than 10 will likely produce shapes
        %that are less accurate. Default value set to 5, which means use every
        %5 points on the grain perimeter in the calculation. Smallest
        %possible increment is 1.
    %meanphif: Specifies what mean intergranular friction angle to use
        % in critical shear stress calculations. In units of degrees.
        %Cannot be more than 89. 
    %stdphif: Specifies what standard deviation of intergranular friction
        %angle to use in critical shear stress calculations. In units of
        %degrees.
    %numphif: Specifies how many intergranular friction
        %angles to sample from a normal distribution.  A higher value (e.g. 2000) 
        %will result in a more accurate distribution that is less subject to
        %random choice. If the distribution will include thetaf of 90 degrees 
        %or larger, then a larger sample size may be needed because many 
        %these values will be automatically deleted from the calculations. 
        %In these cases the actual mean and standard deviation will not 
        %always equal the set values because of truncation at the upper end
        %of the distribution
    %rhos: Assumed or measured average sediment density. In units of kg/m^3
    %porosity of bed sediment. Varies from 0-1. Typical values are between
        %0.4-0.6
    %Cd: Drag coefficient for sediment. 
    %Cl: Lift coefficient for sediment.
    %Cv: Volume correction factor for resisting force calculations to account 
        %for more sediment that is resisting grain motion than just its own buried volume.
        %This should not be adjusted unless you have data to support a different
        %value. See Yager et al. (2018, JGR) for details on this constant.
    %K: von Karmen's constant    
    %rho: Water density in units of kg/m^3
    %g: Gravitational acceleration in units of m^2/s
    %phistep: Interval in phi units to bin grain sizes and resulting
        %protrusion and critical shear stress values. Recommend using half-phi
        %intervals (phistep=-0.5) but whole phi intervals phistep=-1) could be more
        %appropriate if the sample size is low, which will cause few grains in
        %each bin and likely large uncertainties in estimated protrusion and
        %critical shear stress distributions for each grain size bin. phistep
        %values must be negative to keep grain size bins in the gravel range.

        
    % If the .csv file does not exist, use default values for each input variable
 else
    disp('--- USING DEFAULT PRO+ INPUT VARIABLES FROM defineinputs.m');
    %see input variable defintions above
    proinputs.saveoutputs=0;
    proinputs.usepivot=0; 
    proinputs.whichinput=1;
    proinputs.whichks=1; 
    proinputs.whichradius=1; 
    proinputs.setks=0.014;
    proinputs.setradius=0.014;
    proinputs.perimpoints=5;
    proinputs.meanphif=60;
    proinputs.stdphif=15;
    proinputs.numphif=2000;
    proinputs.phistep=-0.5;   
    proinputs.rhos=2650;
    proinputs.rho=1000;
    proinputs.porosity=0.4;
    proinputs.Cd=0.4;
    proinputs.Cl=0.2;
    proinputs.Cv=4.827;
    proinputs.K=0.4;
    proinputs.g=9.81;

end    
%% Folder to save outputs, will overwrite files to avoid saving many files with different names!
if proinputs.saveoutputs==1  
    test = exist('ProPlus/','dir');  
    if test==0 
        mkdir('ProPlus/');  
    end 
end
end
