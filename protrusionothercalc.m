function [medianp,percentile10p]=protrusionothercalc(proinputs,bedstd,b84,ptCloud,grain) 
%% Function to calculate protrusion for each grain on the bed using outputs from G3Point and a Point Cloud 
%Two variables are output from this function: percentile10p, which is the protrusion
%calculated using the 10th percentile of the surrounding bed elevation and
%is used in driving force equations for critical shear stress calculations;
%medianp,which is the protrusion calculated using the 10th percentile of 
%the surrounding bed elevation and is used in resisting force equations 
%for critical shear stress calculations. This code could be modified to
%calcualte protrusion using other inputs such as grain perimeters from
%machine learning applications
%% Determine what search radius defintion to use based on user input
if proinputs.whichradius==1 
    search_radius=bedstd;
elseif proinputs.whichradius==2
    search_radius=b84;
elseif proinputs.whichradius==3
    search_radius=proinputs.setradius;
else
    disp('error, need to enter search radius choice between 1-3')
end
%% CALCULATE PROTRUSION FOR EACH GRAIN
disp('--- CALCULATING PROTRUSION')
%preallocate the output variables
numgrains=size(grain,1);
percentile10p=NaN.*ones(numgrains,1);medianp=NaN.*ones(numgrains,1);

%loop through all identified grains
countgrains=1;counter=1;
    for j=1:numgrains            
        %find the point cloud for each grain using the specified perimeter
        [indin,indon]=inpolygon(ptCloud.Location(:,1),ptCloud.Location(:,2),grain(j).perim(:,1),grain(j).perim(:,2));
        %identify the top of each grain
        graintop=max(ptCloud.Location(indin|indon,3));
        %find the point cloud that does not include the grain of interest:
        %two different point clouds are used, one with 0 elevation for z for radius
        %search, the other with actual elevations for protrusion
        %calculation
        if isempty(graintop)==0
            ptCloud_outside=pointCloud([ptCloud.Location(~indin&~indon,1),ptCloud.Location(~indin&~indon,2),ptCloud.Location(~indin&~indon,3)]);
            ptCloud_search=pointCloud([ptCloud.Location(~indin&~indon,1),ptCloud.Location(~indin&~indon,2),zeros(length(ptCloud_outside.Location(:,3)),1)]);
        
        %search the point cloud to find points that are within a certain radius (x and y) of
        %every proinputs.perimpoint spacing of points on the perimeter of the grain of interest
            indcircle=NaN.*ones(500000,size(grain(j).perim,1));
            for pp=1:proinputs.perimpoints:size(grain(j).perim,1)      
                [indcircled,~] = findNeighborsInRadius(ptCloud_search,[grain(j).perim(pp,1),grain(j).perim(pp,2),0],search_radius);
                indcircle(1:numel(indcircled),pp)=indcircled;
            end
            indcircle=reshape(indcircle,size(indcircle,1)*size(indcircle,2),1);indcircle(isnan(indcircle))=[];
        
        %use these points within the specified radius to calculate
        %protrusion and upstream and downstream protrusion
            ptCloudprotrusion = select(ptCloud_outside,unique(indcircle));
            medianp(countgrains)=graintop-median(ptCloudprotrusion.Location(:,3));
            percentile10p(countgrains)=graintop-prctile(ptCloudprotrusion.Location(:,3),10);
            countgrains=countgrains+1; percom=round((j/numgrains)*100);
            if (percom==10 && counter==1) || (percom==25 && counter==2)...
                || (percom==50 && counter==3) || (percom==75 && counter==4)...
                || (percom==90 && counter==5)...
                ;counter=counter+1;
                fprintf('Protrusion calculation is %d percent complete\n',percom);
            end
            clear indin indon indcircled ptCloudprotrusion
        end
        clear ptCloud_outside ptCloud_search graintop boundind indcircle
    end
end