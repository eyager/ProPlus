function [FrWt_dist,taucr_dist,taustcr_dist,phif_mean_real,phif_std_real,b_dist,pD_dist,pR_dist] = taucrcalc(proinputs,b,b50,b84,bedstd,pD,pR)
%% USER NOTE
%Code to calculate the resisting forces (Yager et al., 2018 JGR), driving forces and
%resulting critical shear stresses for a distribution of grains on the bed.
%Use this code at your own risk and modifications may be needed for your
%specific application. If you encounter errors in the code please tell
%Elowyn Yager (eyager@uidaho.edu).  Code written by Elowyn Yager and last
%modified on 11/15/2023.
%% OUTPUT VARIABLES
%thetaf_mean_real: actual mean of the thetaf distribution employed, which
    %may differ from the originally input thetaf_mean value (see below for
    %details). In units of degrees.
%thetaf_std_real: actual standard deviation of the thetaf distribution employed, which
    %may differ from the originally input thetaf_std value (see below for
    %details). In units of degrees.
%FrWt_dist: estimated weight normalized distribution of estimated submerged resisting
    %forces. This can directly compared with normalized measured
    %resisting forces even if those measured forces were for an unsubmerged
    %bed. Measured resisting forces on unsubmerged beds should be normalized by 
    %the unsubmerged weight. A comparision of measured and estimated normalized
    %resisting forces allows you to adjust the assumed thetaf distribution
    %values until the measured and estimated normalized resisting forces match.
    %This essentially allows for indirect calibration of thetaf in the critical shear
    %stress calculations. Dimensionless.
%taucr_dist: estimated critical shear stress distribution in units of Pa.
%taustcr_dist: estimated critical Shields stress distribution.
    %Dimensionless
%b_dist: grain sizes that correspond to each estimated critical shear
    %stress. In units of meters
%pD_dist: driving protrusion that corresponds to each estimated critical
    %shear stress, this can be larger than the input pD distribution because
    %assumed thetaf and thetap distributions cause the critical shear stress
    %distribution to have more values. This pD_dist is mostly for
    %reference.In units of meters.
%pR_dist: resisting protrusion that corresponds to each estimated critical 
    %shear stress, this can be larger than the input pR distribution because
    %assumed thetaf and thetap distributions cause the critical shear stress
    %distribution to have more values. This pR_dist is mostly for reference
    %but some pR values were changed in the code below because of requirements
    %in the critical shear stress calculations. This may cause the
    %input pR and pR_dist distributions to differ but pR_dist is provided 
    %mostly for reference.In units of meters.

%% CALCULATE BASIC GRAIN VARIABLES
disp('--- CALCULATING FORCES AND CRITICAL SHEAR STRESSES');
r=b./2; %grain radius (r) is half particle diameter
Vgrain=(4./3).*pi.*(r.^3); %volume of each particle (assumed spheres)
Wt=(proinputs.rhos-proinputs.rho).*proinputs.g.*Vgrain; %submerged particle weight
indpR=find(pR>b);pR(indpR)=b(indpR);%set the resisting protrusion equal to the grain diameter 
 %if protrusion exceeds the grain diameter 
pR(pR<=0)=0;%set resisting protrusion equalto zero if protrusion is negative.
%% CALCULATE BASIC BED/FLOW VARIABLES
pore=1-proinputs.porosity; %One minus the bed porosity (packing fraction)

%calculate the thetaf distribution in radians using the user defined mean, 
%standard deviation,and number of points to sample from normal distribution
phif=normrnd(pi.*proinputs.meanphif./180,pi.*proinputs.stdphif./180,proinputs.numphif,1);
phif(tan(phif)<0)=NaN; %eliminate any thetaf >= 90 degrees
phif_mean_real=nanmean(180.*phif./pi);%calculate the actual mean thetaf after removing values
phif_std_real=nanstd(180.*phif./pi);%calculate the actual std thetaf after removing values

%set the percentiles (perc) of the pivot angle (thetap) distribution that
%will be sampled from the Kirchner et al. 1990 equation (see Yager et al.
%2018 JGR for details) to obtain the calculated resisting force.
%Recommended that percentiles steps of at least 1 are used to obtain a large enough
%sample size and that the full range of possible percentiles are employed
%(1-100). pivot angles are not used if proinputs.usepivot=0
if proinputs.usepivot==0
    perc=NaN;
elseif proinputs.usepivot==1
    perc=1:1:100;
else
    disp('error, need to enter usepivot choice of 0 or 1')
end

%determine which roughness length (ks) defintion to use based on user input
if proinputs.whichks==1
    ks=bedstd; %use standard deviation of bed elevations
elseif proinputs.whichks==2
    ks=b84; %use 84th percentile of GSD
elseif proinputs.whichks==3
    ks=proinputs.setks;    %use user set value
else
    disp('error, need to enter ks choice between 1-3')
end
%% PREDEFINE VARIABLES FOR LOOP
%many of these variables will not be exported but are saved here to avoid
%variable growing in loop,in case any troubleshooting needs to occur, or if
%these intermediate variables might want to be exported by a user

fpD=NaN.*ones(length(b),1);%lower bound of velocity profile for lift
fp=NaN.*ones(length(b),1);%upper bound of velocity profile for lift
intfx=NaN.*ones(length(b),1);%integral in Lamb et al. (2017) velocity profile
phip=NaN.*ones(length(perc),length(b));%pivot angle distribution
Vb=NaN.*ones(length(b),1);%buried grain volume
Vo=NaN.*ones(length(b),1);%overlying sediment volume
Fs=NaN.*ones(length(perc),length(b));%force due to burial
Fg=NaN.*ones(length(perc),length(b));%force due tograin weight
Fsed=NaN.*ones(length(phif),length(b));%force due to intergranular friction
Fr=NaN.*ones(length(perc),length(b),length(phif));%total resisting force
FrWt=NaN.*ones(length(perc),length(b),length(phif));%normalized total resisting force
taucr=NaN.*ones(length(perc),length(b),length(phif));%critical shear stress (in Pa)
taustcr=NaN.*ones(length(perc),length(b),length(phif));%critical Shields stress
b_all=NaN.*ones(length(perc),length(b),length(phif));%grain sizes that correspond to each calculated critical shear stress
pD_all=NaN.*ones(length(perc),length(b),length(phif));%driving protrusion that corresponds to each calculated critical shear stress
pR_all=NaN.*ones(length(perc),length(b),length(phif));%resisting protrusion that corresponds to each calculated critical shear stress

%% CALCULATE RESISTING AND DRIVING FORCES AND TAUCR* FOR EACH GRAIN
for i=1:length(b) %first for loop for calculations that only rely on each grain (not thetap or thetaf)
    
    %if statement to determine if half or more of the grain is buried (pR<r), which is
    %needed to determine the overlying sediment volume for grain resistance calculation (Vo; eqn 4 in Yager
    %et al., 2018; note that there were typos in the published equation--the equation in this code
    %is correct and is what all results in that paper were calculated using). If half is not buried, Vo is zero.
    if pR(i)<r(i)
        Vo(i)=(((b(i).^3)./2)-((b(i).^2).*pR(i)))-(pi.*(((b(i).^3)./12)+((pR(i).^3)./3)-(((pR(i).^2).*b(i))./2)));
    else 
        Vo(i)=0;
    end
    
    %calculate the buried grain volume (Vb; eqn 6 in Yager et al. 2018)
    Vb(i)=(pi.*((b(i)-pR(i)).^2).*((b(i)./2)-((b(i)-pR(i))./3)));
    
    %if statements to determine f(z) (velocity profile function) on the top (fp) and 
    %bottom (fpD) of the grain for lift force calculations. See eqn 11 in
    %Lamb et al. 2017 for details. u(z)=0 at z=0, where z is the distance above the bed, 
    %and where z=0 at the 10th percentile of the surrounding bed elevation.
    %pD-b (grain bottom) or pD (grain top) are substituted into eqn 11 for
    %z. Any grain without full protrusion will have a velocity of 0 at the
    %grain bottom because the grain bottom is below the surrounding bed
    %elevation. Any grain with negative protrusion will have 0 velocity on
    %the entirety of the particle because the particle is fully below the
    %surrounding bed elevation.
    if (pD(i)-b(i))<0 || pD(i)<0 
        fpD(i)=0;
    else
        fpD(i)=log(1+((30.*(pD(i)-b(i)))./ks));
    end
    
    if pD(i)<0
        fp(i)=0;
    else
        fp(i)=log(1+((30.*(pD(i)))./ks));
    end

    %set up numerical integration function (fun1) for the velocity profile 
    %integral in eqn 12 of Kirchner et al. (1990) 
    %but where the velcoity profile part of the integral is replaced using
    %the integral part of eqn 11 in Lamb et al. 2017. The first part of the
    %fun1 accounts for grain area exposed to the flow and the second part
    %(log term) accounts for the velocity profile
     fun1=@(z) (((b(i).^2)-(((2.*z)-((2.*pD(i))-(b(i)))).^2)).^0.5).*((log(1+((30.*z)./ks))).^2); 
     
    %conditional statments to determine the correct limits of integration depending
    %on the particle protrusion.  When pD is less than
    %0, the particle is not exposed to any flow and will not have a
    %solution for the integral (intfx). The upper limit of integration is 
    %pD. If pD is less than b, 
    %the bottom of the particle sits below the surrounding bed elevation
    %and therefore the lower limit of integration is at z=0. If pD is greater
    %than b, then the bottom of the particle sits above the surrounding bed 
    %elevation adn lower limit of integration is pD-b. 
    if (pD(i)-b(i))<0 && pD(i)>0
          intfx(i)=integral(@(z)fun1(z),0,pD(i));
    elseif pD(i)<0
        intfx(i)=NaN;
    else 
        intfx(i)=integral(@(z)fun1(z),pD(i)-b(i),pD(i));
    end
    clear fun1
    
    %second loop to calculate variables that depend on thetap; variables
    %are calculated for every combination of thetap for every grain except
    %when proinputs.usepivot=0, then pivot angle is not included
    for j=1:length(perc)
        if isnan(perc(j))
        phip(j,i)=0.785398163; % results in tan(thetap)=1 and pivot angle is effectively 
        %removed from the calculations
        else
        phip(j,i)=(pi.*((30+0.5.*perc(j)).*((b(i)./b50).^(-0.3)))./180);%pivot angle
        %calculated using equation in Kirchner et al., 1990. Given in
        %radians.
            if tan(phip(j,i))<0 
                phip(j,i)=NaN; %remove any thetap >= 90 degrees
            end
        end
       
        %calculate the submerged force caused by the sediment partly burying the grain of interest
        %(Fs, eqn 3 in Yager et al., 2018). Note that in this and all
        %subsequent force equations, the effects of grain submergence under
        %water are include (rhos-rho) instead of assuming grains on a dry
        %bed, which is in the equations of Yager et al.  If dry bed
        %calculations are needed for direct comparison to measured values, then
        %the '-rho' needs to be deleted from all grain resistance
        %equations but this is not appropriate for taucr calculations.   
        Fs(j,i)=proinputs.g.*(proinputs.rhos-proinputs.rho).*Vo(i).*pore.*tan(phip(j,i));  
        
        %the force from the sumberged grain weight (eqn 1 in Yager et al., 2018) 
        %modified by the pivot angle. See notes above about use of rhos-rho above.
        Fg(j,i)=Wt(i).*tan(phip(j,i));
    
        %third loop to calculate variables that depend on thetaf
        for k=1:length(phif)
            phif(k)=normrnd(phimean,pi.*proinputs.stdphif./180,proinputs.numphif,1);
            phif(tan(phif(k))<0)=NaN; %eliminate any thetaf >= 90 degrees  
            
            
            % submerged force from intergranular friction of surrounding sediment (eqn 5 in
            %Yager et al., 2018; named Fd in paper but has been renamed as
            %Fsed here to avoid confusion with the drag force below).
            Fsed(k,i)=proinputs.Cv.*proinputs.g.*(proinputs.rhos-proinputs.rho).*Vb(i).*pore.*tan(phif(k));
            
            %Submerged total resisting force is calculated for every combination 
            %of thetap, thetaf, and each grain.
            Fr(j,i,k)=(Fg(j,i)+Fsed(k,i)+Fs(j,i));
            
            %Total submerged resisting force normalized by the submerged grain weight;
            FrWt(j,i,k)=Fr(j,i,k)./Wt(i); 
            
            %critical shear stress (taucr; in Pa) from eqn 12 of Kirchner et al.,
            %1990 with modficiation by Yager et al., 2018 to include total Fr
            %rather than just grain weight and pivot angle
            taucr(j,i,k)=(Fr(j,i,k)./(((proinputs.Cd.*intfx(i))./(2.*(proinputs.K.^2)))+(tan(phip(j,i))...
                .*(pi./8).*(proinputs.Cl/(proinputs.K.^2)).*(b(i).^2).*((fp(i).^2)-(fpD(i)).^2))));   
            
           %critical Shields stress
             taustcr(j,i,k)=taucr(j,i,k)./((proinputs.rhos-proinputs.rho).*proinputs.g.*b(i));
             
           %grain size and protrusion values that correspond to each
           %critical shear stress
             b_all(j,i,k)=b(i);pD_all(j,i,k)=pD(i);pR_all(j,i,k)=pR(i);
        end
    end
end
clear i j k

%% RESHAPE OUTPUTS
%entire distribution of values for the entire sample combining all grain sizes, pivot angles, and
%intergranular friction angles.  This is what is output from the function.
FrWt_dist=reshape(FrWt,size(FrWt,1).*size(FrWt,2).*size(FrWt,3),1);%estimated normalized submerged resisting force distribution
taucr_dist=reshape(taucr,size(taucr,1).*size(taucr,2).*size(taucr,3),1);%estimated critical shear stress distribution
taustcr_dist=reshape(taustcr,size(taustcr,1).*size(taustcr,2).*size(taustcr,3),1);%estimated critical Shields stress distribution
b_dist=reshape(b_all,size(b_all,1).*size(b_all,2).*size(b_all,3),1);%b that correspond to each estimated critical Shields stress 
pD_dist=reshape(pD_all,size(pD_all,1).*size(pD_all,2).*size(pD_all,3),1);%pD that correspond to each estimated critical Shields stress 
pR_dist=reshape(pR_all,size(pR_all,1).*size(pR_all,2).*size(pR_all,3),1);%pR that correspond to each estimated critical Shields stress 
end








