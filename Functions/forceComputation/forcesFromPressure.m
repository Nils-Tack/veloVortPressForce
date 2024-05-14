function [summaryData,allLocalForceData] = forcesFromPressure(D,path,opts,plottingOptions,firstLast)
% various options to change appearance of quiver plots
arrowThickness = 1;                 % thickness of quiver arrows 
arrowHead = 0.15;                   % arrow head proportion relative to length of quiver
pullThrustCol = [0.45 0.68 0.8];    % color code of pull thrust force vectors
pushThrustCol = [0.8 0.5 0.5];      % color code of push thrust force vectors
pullDragCol = [0 0.2 0.5];          % color code of pull drag force vectors
pushDradCol = [0.5 0 0.1];          % color code of push drag force vectors

% parameters to control the scalebar and reference vector
scalebarLength = 0.02; % scale bar in meter (2cm), change as needed
arrowLength = 0.008; % reference scale arrow = 0.08 mN cm-1 = 0.8 mN m

% Set temporal sale 
deltat = opts.increment*opts.pressureDeltaT; % time between frames for force calculation

% Load one image for scaling
Itemp = importdata(quickfilepath(D.images(1)));
I = Itemp(:,:,1); % select only one of the three layers of the image;

% Load heading data
if opts.heading == 1 % fixed heading
heading = opts.fixedHeading;
elseif opts.heading == 2 % changing heading
heading = importdata(quickfilepath(D.heading(1)));
end

% Set figure and force data export path
if firstLast(1)~=firstLast(2)
PathLocalForces  = fullfile(path.forces,'local force-power');
mkdir(PathLocalForces);
end

if opts.exportFigures == 1 && firstLast(1)~=firstLast(2)
PathForcesFigures = fullfile(path.forces,'figures'); 
mkdir(PathForcesFigures) % make new directory for force figures
end

% Figure options for testing vs. all frames
if firstLast(1)~=firstLast(2)
   if strcmp('all',plottingOptions)
        figure('units','centimeters','Position',[1 1 30 10]);
   else
        figure('units','centimeters','Position',[1 1 20 20]);
   end
end

advanceFr = 1;                              % Initiate variable to advance  index of force data to be saved
summaryData = zeros(firstLast(2),13);       % Initiate summary data table to be output
allLocalForceData = cell(firstLast(2),1);   % Initiate cell array to store local force and power data


for fr = firstLast(1):opts.increment:firstLast(2)
% Read files
    % Image
    if opts.useRawImage == 1
        if fr == firstLast(1)
            Itemp = importdata(quickfilepath(D.images(fr)));
            I = Itemp(:,:,1); % select only one of the three layers of the image;
            % I = D_tif(:,:,fr); % Import image for plotting purposes
        else
            Itemp = importdata(quickfilepath(D.images(fr-1)));
            I = Itemp(:,:,1); % select only one of the three layers of the image;
        end
    end
    
    % Mask file - set to desired number of points along the boundary - uses curvspace
    blank = importdata(quickfilepath(D.outlines(fr))); % import outlines blanking coordinates
    blank = curvspace(blank,opts.npoints); % number of points to re-define outline is wanted
    blank = alignPoints(blank); % always align array with anteriormost point of animal as first element ---> this is where centerlines coule be useful to find the snout and define the begining of the loop making up the outline
    boundx = blank(:,1);  % set x coordinates in m
    boundy = blank(:,2);  % set y coordinates in m
   
    % Pressure file
    press = importdata(quickfilepath(D.pressure(fr))); % import outlines blanking coordinates
    [Xpress,Ypress] = meshgrid(unique(press(:,1),'stable'),unique(press(:,2),'stable'));
    PP = reshape(press(:,7),[size(Xpress,1),size(Xpress,2)]);

% Test plot
if firstLast(1)==firstLast(2)
figure('units','centimeters','Position',[1 1 30 15]);
subplot(1,2,1);
hold on
if opts.useRawImage == 1
imagesc([0,(size(I,2)-1)/opts.scale],[0,(size(I,1)-1)/opts.scale],flipud(I)); colormap('gray') % plot image --------------------> flipud I

end
plot(boundx,boundy,'-r','Linewidth',1.5);
axis equal
title(sprintf('Outline - Frame %i',fr))

ax(1) = subplot(1,2,2);
axisMagnitude = 10; % pressure magnitude
hold on
contourf(Xpress,Ypress,PP,150,'edgecolor','none') % plot pressure fields
clim([-axisMagnitude axisMagnitude]) % pressure field min and max magnitude
colormap(ax(1),customcolormap_preset('orange-white-purple'))    
patch(boundx,boundy,'k') % overlay outlines
hcb = colorbar;
hcb.Label.String = 'relative pressure (Pa)';
axis equal
title(sprintf('Pressure field - Frame %i',fr))
end

% Calculate the coordinates of the center of each segment where the forces should be calculated/plotted
boundx2 = mean([boundx circshift(boundx,1)],2);
boundy2 = mean([boundy circshift(boundy,1)],2);

% calculate pressure at body boundary
surfpress = griddata(Xpress,Ypress,double(PP),boundx2,boundy2);

% Calculate length of each surface segment (equivalent to surface area per unit depth) between surface points. 
surfdx = boundx - circshift(boundx,1);   % compute x-distance between surface points
surfdy = boundy - circshift(boundy,1);   % compute y-distance between surface points
darea = sqrt(surfdx.^2 + surfdy.^2);     % compute surface area (per unit depth) between surface points; also equivalent to segment length

% Compute normal unit vectors (assuming that the animal is swimming DOWN)
surfunitnormx = -surfdy./darea; % compute x-component of vector normal to surface - lateral component
surfunitnormy = surfdx./darea; % compute y-component of vector normal to surface - axial component

% Perform decomposition to compute forces for the origin/axes rotated CW relative to the +y image axis by theta
% angle in the direction of swimming
% WARNING: not the data to use for plotting purposes. It only allows to
% compute the forces for a specific orientation

if opts.heading == 1 % fixed heading
theta = opts.fixedHeading;
elseif opts.heading == 2 % changing heading
theta = heading(fr);
end

surfunitnormx2 = surfunitnormx*cosd(theta)-surfunitnormy*sind(theta); % lateral component
surfunitnormy2 = surfunitnormx*sind(theta)+surfunitnormy*cosd(theta); % axial component

% thrust components for orientation of the vectors when plotting
Tx = surfunitnormy2*sind(theta);
Ty = surfunitnormy2*cosd(theta);

% Lateral components for orientation of the vectors when plotting
Lx = surfunitnormx2*cosd(theta);
Ly = -surfunitnormx2*sind(theta);


% Identify surface points where pressure is pushing/pulling in the same/oposite direction of swimming
indpospull = find(surfunitnormy2 < 0 & surfpress < 0);   % find surface points where low pressure is pulling animal forward
indpospush = find(surfunitnormy2 > 0 & surfpress > 0);   % find surface points where high pressure is pushing animal forward
indnegpull = find(surfunitnormy2 > 0 & surfpress < 0);   % find surface points where low pressure is pulling animal backward
indnegpush = find(surfunitnormy2 < 0 & surfpress > 0);   % find surface points where high pressure is pushing animal backward

% Identify surface points where pressure is pushing/pulling laterally (lateral component x)
indleftpull = find(surfunitnormx2 < 0 & surfpress < 0);   % find surface points where low pressure is pulling animal forward
indleftpush = find(surfunitnormx2 < 0 & surfpress > 0);   % find surface points where high pressure is pushing animal forward
indrightpull = find(surfunitnormx2 > 0 & surfpress < 0);   % find surface points where low pressure is pulling animal backward
indrightpush = find(surfunitnormx2 > 0 & surfpress > 0);   % find surface points where high pressure is pushing animal backward

%% Plot vectors for pressure according to the 4 categories above
if strcmp('axial',plottingOptions)
    if firstLast(1)==firstLast(2)
    figure('units','centimeters','Position',[1 1 20 20]);
    end
    
    hold on
        if opts.useRawImage == 0 
           patch(boundx,boundy,'k')
        else
           imagesc([0,(size(I,2)-1)/opts.scale],[0,(size(I,1)-1)/opts.scale],flipud(I)); colormap('gray') % plot image
        end

    % Quiver plot of forces along body
    % WARNING: use of 'opts.forceScale' for plotting purposes only
    quiver(boundx2(indpospull,:),boundy2(indpospull,:),surfpress(indpospull,:).*Tx(indpospull,:).*darea(indpospull,:)*opts.forceScale,surfpress(indpospull,:).*Ty(indpospull,:).*darea(indpospull,:)*opts.forceScale,0,'Color',pullThrustCol,'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);   % original color was [0 1 1]
    quiver(boundx2(indpospush,:),boundy2(indpospush,:),surfpress(indpospush,:).*Tx(indpospush,:).*darea(indpospush,:)*opts.forceScale,surfpress(indpospush,:).*Ty(indpospush,:).*darea(indpospush,:)*opts.forceScale,0,'Color',pushThrustCol,'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);   % original color was [1 0.5 1]
    quiver(boundx2(indnegpull,:),boundy2(indnegpull,:),surfpress(indnegpull,:).*Tx(indnegpull,:).*darea(indnegpull,:)*opts.forceScale,surfpress(indnegpull,:).*Ty(indnegpull,:).*darea(indnegpull,:)*opts.forceScale,0,'Color',pullDragCol,'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);   % original color was [0 0 1]
    quiver(boundx2(indnegpush,:),boundy2(indnegpush,:),surfpress(indnegpush,:).*Tx(indnegpush,:).*darea(indnegpush,:)*opts.forceScale,surfpress(indnegpush,:).*Ty(indnegpush,:).*darea(indnegpush,:)*opts.forceScale,0,'Color',pushDradCol,'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);   % original color was [1 0 0]
    
        if opts.useRawImage == 0 
            %---Scalebar---
            % if the axes are turned on it may not be necessary. But if ensures the
            % figure is properly formatted all the time
            % scalebarLength = 0.02; % scale (2cm), change as needed
            scaleBarXmax = ((size(I,2)-1)/opts.scale)*0.97; % extremity at 97% of the right edge of the frame
            scaleBarXmin = scaleBarXmax-scalebarLength; %  
            scaleBarY = ((size(I,1)-1)/opts.scale)*0.03; % extremity at 3% of the bottom edge of the frame
            plot([scaleBarXmin; scaleBarXmax], [scaleBarY; scaleBarY], '-k', 'LineWidth', 3) % change the color of the bar as needed. '-k' for black, '-w' for white
            text(scaleBarXmin,(size(I,1)-1)/opts.scale*0.05,sprintf('%i cm',scalebarLength*100),'Color',[0 0 0]) % change the text and its color as needed
            
            %---Reference vector---
            p1 = [((size(I,2)-1)/opts.scale)*0.03 ((size(I,1)-1)/opts.scale)*0.03]; % First Point at 3% of max y axis and x axis (x and y of origin of scale vector)
            quiver(p1(1),p1(2),0,opts.forceScale*arrowLength,0,'MaxHeadSize',0.7,'LineWidth', 2,'Color', 'k') % reference vector scale
            text(((size(I,2)-1)/opts.scale)*0.035,((size(I,1)-1)/opts.scale)*0.04,sprintf('%.3f mN cm^{-1}',arrowLength),'Color',[0 0 0])
        else
            %---Scalebar---
            % if the axes are turned on it may not be necessary. But if ensures the
            % figure is properly formatted all the time
            % scalebarLength = 0.02; % scale (2cm), change as needed
            scaleBarXmax = ((size(I,2)-1)/opts.scale)*0.97; % extremity at 97% of the right edge of the frame
            scaleBarXmin = scaleBarXmax-scalebarLength; %  
            scaleBarY = ((size(I,1)-1)/opts.scale)*0.03; % extremity at 3% of the bottom edge of the frame
            plot([scaleBarXmin; scaleBarXmax], [scaleBarY; scaleBarY], '-r', 'LineWidth', 3) % change the color of the bar as needed. '-k' for black, '-w' for white
            text(scaleBarXmin,(size(I,1)-1)/opts.scale*0.05,sprintf('%i cm',scalebarLength*100),'Color',[1 0 0]) % change the text and its color as needed

            %---Reference vector---
            p1 = [((size(I,2)-1)/opts.scale)*0.03 ((size(I,1)-1)/opts.scale)*0.03]; % First Point at 3% of max y axis and x axis (x and y of origin of scale vector)
            quiver(p1(1),p1(2),0,opts.forceScale*arrowLength,0,'MaxHeadSize',0.7,'LineWidth', 2,'Color', 'r') % reference vector scale
            text(((size(I,2)-1)/opts.scale)*0.035,((size(I,1)-1)/opts.scale)*0.04,sprintf('%.3f mN cm^{-1}',arrowLength),'Color',[1 0 0])
        end

    title(sprintf('Axial force - frame %i',fr))
    axis([0,(size(I,2)-1)/opts.scale, 0,(size(I,1)-1)/opts.scale])
    axis off
    axis equal
    formatFigure
    
elseif strcmp('lateral',plottingOptions)
        if firstLast(1)==firstLast(2)
        figure('units','centimeters','Position',[1 1 20 20]);
        end
        
        hold on
            if opts.useRawImage == 0 
               patch(boundx,boundy,'k')
            else
               imagesc([0,(size(I,2)-1)/opts.scale],[0,(size(I,1)-1)/opts.scale],flipud(I)); colormap('gray') % plot image
            end

    quiver(boundx2(indleftpull,:),boundy2(indleftpull,:),surfpress(indleftpull,:).*Lx(indleftpull,:).*darea(indleftpull,:)*opts.forceScale,surfpress(indleftpull,:).*Ly(indleftpull,:).*darea(indleftpull,:)*opts.forceScale,0,'Color','#562689','LineWidth',arrowThickness,'MaxHeadSize',arrowHead);   % original color was [0 1 1]
    quiver(boundx2(indleftpush,:),boundy2(indleftpush,:),surfpress(indleftpush,:).*Lx(indleftpush,:).*darea(indleftpush,:)*opts.forceScale,surfpress(indleftpush,:).*Ly(indleftpush,:).*darea(indleftpush,:)*opts.forceScale,0,'Color','#b35807','LineWidth',arrowThickness,'MaxHeadSize',arrowHead);
    quiver(boundx2(indrightpull,:),boundy2(indrightpull,:),surfpress(indrightpull,:).*Lx(indrightpull,:).*darea(indrightpull,:)*opts.forceScale,surfpress(indrightpull,:).*Ly(indrightpull,:).*darea(indrightpull,:)*opts.forceScale,0,'Color','#562689','LineWidth',arrowThickness,'MaxHeadSize',arrowHead);
    quiver(boundx2(indrightpush,:),boundy2(indrightpush,:),surfpress(indrightpush,:).*Lx(indrightpush,:).*darea(indrightpush,:)*opts.forceScale,surfpress(indrightpush,:).*Ly(indrightpush,:).*darea(indrightpush,:)*opts.forceScale,0,'Color','#b35807','LineWidth',arrowThickness,'MaxHeadSize',arrowHead);

        if opts.useRawImage == 0 
            %---Scalebar---
            % if the axes are turned on it may not be necessary. But if ensures the
            % figure is properly formatted all the time
            % scalebarLength = 0.02; % scale (2cm), change as needed
            scaleBarXmax = ((size(I,2)-1)/opts.scale)*0.97; % extremity at 97% of the right edge of the frame
            scaleBarXmin = scaleBarXmax-scalebarLength; %  
            scaleBarY = ((size(I,1)-1)/opts.scale)*0.03; % extremity at 3% of the bottom edge of the frame
            plot([scaleBarXmin; scaleBarXmax], [scaleBarY; scaleBarY], '-k', 'LineWidth', 3) % change the color of the bar as needed. '-k' for black, '-w' for white
            text(scaleBarXmin,(size(I,1)-1)/opts.scale*0.05,sprintf('%i cm',scalebarLength*100),'Color',[0 0 0]) % change the text and its color as needed

            %---Reference vector---
            p1 = [((size(I,2)-1)/opts.scale)*0.03 ((size(I,1)-1)/opts.scale)*0.03]; % First Point at 3% of max y axis and x axis (x and y of origin of scale vector)
            quiver(p1(1),p1(2),0,opts.forceScale*arrowLength,0,'MaxHeadSize',0.7,'LineWidth', 2,'Color', 'k') % reference vector scale
            text(((size(I,2)-1)/opts.scale)*0.035,((size(I,1)-1)/opts.scale)*0.04,sprintf('%.3f mN cm^{-1}',arrowLength),'Color',[0 0 0])
        else
            %---Scalebar---
            % if the axes are turned on it may not be necessary. But if ensures the
            % figure is properly formatted all the time
            % scalebarLength = 0.02; % scale (2cm), change as needed
            scaleBarXmax = ((size(I,2)-1)/opts.scale)*0.97; % extremity at 97% of the right edge of the frame
            scaleBarXmin = scaleBarXmax-scalebarLength; %  
            scaleBarY = ((size(I,1)-1)/opts.scale)*0.03; % extremity at 3% of the bottom edge of the frame
            plot([scaleBarXmin; scaleBarXmax], [scaleBarY; scaleBarY], '-w', 'LineWidth', 3) % change the color of the bar as needed. '-k' for black, '-w' for white
            text(scaleBarXmin,(size(I,1)-1)/opts.scale*0.05,sprintf('%i cm',scalebarLength*100),'Color',[1 1 1]) % change the text and its color as needed

            %---Reference vector---
            p1 = [((size(I,2)-1)/opts.scale)*0.03 ((size(I,1)-1)/opts.scale)*0.03]; % First Point at 3% of max y axis and x axis (x and y of origin of scale vector)
            quiver(p1(1),p1(2),0,opts.forceScale*arrowLength,0,'MaxHeadSize',0.7,'LineWidth', 2,'Color', 'w') % reference vector scale
            text(((size(I,2)-1)/opts.scale)*0.035,((size(I,1)-1)/opts.scale)*0.04,sprintf('%.3f mN cm^{-1}',arrowLength),'Color',[1 1 1])
        end

    title(sprintf('Lateral force - frame %i',fr))
    axis([0,(size(I,2)-1)/opts.scale, 0,(size(I,1)-1)/opts.scale])
    axis off
    axis equal
    formatFigure
 
elseif strcmp('both',plottingOptions)
    if firstLast(1)==firstLast(2)    
    figure('units','centimeters','Position',[1 1 20 20]);
    end
    
        hold on
            if opts.useRawImage == 0 
               patch(boundx,boundy,'k')
            else
               imagesc([0,(size(I,2)-1)/opts.scale],[0,(size(I,1)-1)/opts.scale],flipud(I)); colormap('gray') % plot image

            end
    % plot(boundx,boundy,'.r')
    quiver(boundx2(indpospull,:),boundy2(indpospull,:),opts.forceScale*surfpress(indpospull,:).*surfunitnormx(indpospull,:).*darea(indpospull,:),opts.forceScale*surfpress(indpospull,:).*surfunitnormy(indpospull,:).*darea(indpospull,:),0,'Color',[0.45 0.68 0.8],'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);
    quiver(boundx2(indpospush,:),boundy2(indpospush,:),opts.forceScale*surfpress(indpospush,:).*surfunitnormx(indpospush,:).*darea(indpospush,:),opts.forceScale*surfpress(indpospush,:).*surfunitnormy(indpospush,:).*darea(indpospush,:),0,'Color',[0.8 0.5 0.5],'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);
    quiver(boundx2(indnegpull,:),boundy2(indnegpull,:),opts.forceScale*surfpress(indnegpull,:).*surfunitnormx(indnegpull,:).*darea(indnegpull,:),opts.forceScale*surfpress(indnegpull,:).*surfunitnormy(indnegpull,:).*darea(indnegpull,:),0,'Color',[0 0.2 0.5],'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);
    quiver(boundx2(indnegpush,:),boundy2(indnegpush,:),opts.forceScale*surfpress(indnegpush,:).*surfunitnormx(indnegpush,:).*darea(indnegpush,:),opts.forceScale*surfpress(indnegpush,:).*surfunitnormy(indnegpush,:).*darea(indnegpush,:),0,'Color',[0.5 0 0.1],'LineWidth',arrowThickness,'MaxHeadSize',arrowHead)

        if opts.useRawImage == 0 
            %---Scalebar---
            % if the axes are turned on it may not be necessary. But if ensures the
            % figure is properly formatted all the time
            % scalebarLength = 0.02; % scale (2cm), change as needed
            scaleBarXmax = ((size(I,2)-1)/opts.scale)*0.97; % extremity at 97% of the right edge of the frame
            scaleBarXmin = scaleBarXmax-scalebarLength; %  
            scaleBarY = ((size(I,1)-1)/opts.scale)*0.03; % extremity at 3% of the bottom edge of the frame
            plot([scaleBarXmin; scaleBarXmax], [scaleBarY; scaleBarY], '-k', 'LineWidth', 3) % change the color of the bar as needed. '-k' for black, '-w' for white
            text(scaleBarXmin,(size(I,1)-1)/opts.scale*0.05,sprintf('%i cm',scalebarLength*100),'Color',[0 0 0]) % change the text and its color as needed

            %---Reference vector---
            p1 = [((size(I,2)-1)/opts.scale)*0.03 ((size(I,1)-1)/opts.scale)*0.03]; % First Point at 3% of max y axis and x axis (x and y of origin of scale vector)
            quiver(p1(1),p1(2),0,opts.forceScale*arrowLength,0,'MaxHeadSize',0.7,'LineWidth', 2,'Color', 'k') % reference vector scale
            text(((size(I,2)-1)/opts.scale)*0.035,((size(I,1)-1)/opts.scale)*0.04,sprintf('%.3f mN cm^{-1}',arrowLength),'Color',[0 0 0])
        else
            %---Scalebar---
            % if the axes are turned on it may not be necessary. But if ensures the
            % figure is properly formatted all the time
            % scalebarLength = 0.02; % scale (2cm), change as needed
            scaleBarXmax = ((size(I,2)-1)/opts.scale)*0.97; % extremity at 97% of the right edge of the frame
            scaleBarXmin = scaleBarXmax-scalebarLength; %  
            scaleBarY = ((size(I,1)-1)/opts.scale)*0.03; % extremity at 3% of the bottom edge of the frame
            plot([scaleBarXmin; scaleBarXmax], [scaleBarY; scaleBarY], '-w', 'LineWidth', 3) % change the color of the bar as needed. '-k' for black, '-w' for white
            text(scaleBarXmin,(size(I,1)-1)/opts.scale*0.05,sprintf('%i cm',scalebarLength*100),'Color',[1 1 1]) % change the text and its color as needed

            %---Reference vector---
            p1 = [((size(I,2)-1)/opts.scale)*0.03 ((size(I,1)-1)/opts.scale)*0.03]; % First Point at 3% of max y axis and x axis (x and y of origin of scale vector)
            quiver(p1(1),p1(2),0,opts.forceScale*arrowLength,0,'MaxHeadSize',0.7,'LineWidth', 2,'Color', 'w') % reference vector scale
            text(((size(I,2)-1)/opts.scale)*0.035,((size(I,1)-1)/opts.scale)*0.04,sprintf('%.3f mN cm^{-1}',arrowLength),'Color',[1 1 1])
        end

        title(sprintf('Force vectors - frame %i',fr))
        axis([0,(size(I,2)-1)/opts.scale, 0,(size(I,1)-1)/opts.scale])
        axis off
        axis equal
        formatFigure
    
elseif strcmp('all',plottingOptions)
    if firstLast(1)==firstLast(2)
    figure('units','centimeters','Position',[1 1 30 10]); 
    end
    
subplot(1,3,1)
    hold on
        if opts.useRawImage == 0 
           patch(boundx,boundy,'k')
        else
           imagesc([0,(size(I,2)-1)/opts.scale],[0,(size(I,1)-1)/opts.scale],flipud(I)); colormap('gray') % plot image
        end

    % Quiver plot of forces along body
    % WARNING: use of 'opts.forceScale' for plotting purposes only
    quiver(boundx2(indpospull,:),boundy2(indpospull,:),surfpress(indpospull,:).*Tx(indpospull,:).*darea(indpospull,:)*opts.forceScale,surfpress(indpospull,:).*Ty(indpospull,:).*darea(indpospull,:)*opts.forceScale,0,'Color',pullThrustCol,'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);   % original color was [0 1 1]
    quiver(boundx2(indpospush,:),boundy2(indpospush,:),surfpress(indpospush,:).*Tx(indpospush,:).*darea(indpospush,:)*opts.forceScale,surfpress(indpospush,:).*Ty(indpospush,:).*darea(indpospush,:)*opts.forceScale,0,'Color',pushThrustCol,'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);   % original color was [1 0.5 1]
    quiver(boundx2(indnegpull,:),boundy2(indnegpull,:),surfpress(indnegpull,:).*Tx(indnegpull,:).*darea(indnegpull,:)*opts.forceScale,surfpress(indnegpull,:).*Ty(indnegpull,:).*darea(indnegpull,:)*opts.forceScale,0,'Color',pullDragCol,'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);   % original color was [0 0 1]
    quiver(boundx2(indnegpush,:),boundy2(indnegpush,:),surfpress(indnegpush,:).*Tx(indnegpush,:).*darea(indnegpush,:)*opts.forceScale,surfpress(indnegpush,:).*Ty(indnegpush,:).*darea(indnegpush,:)*opts.forceScale,0,'Color',pushDradCol,'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);   % original color was [1 0 0]

        if opts.useRawImage == 0 
            %---Scalebar---
            % if the axes are turned on it may not be necessary. But if ensures the
            % figure is properly formatted all the time
            % scalebarLength = 0.02; % scale (2cm), change as needed
            scaleBarXmax = ((size(I,2)-1)/opts.scale)*0.97; % extremity at 97% of the right edge of the frame
            scaleBarXmin = scaleBarXmax-scalebarLength; %  
            scaleBarY = ((size(I,1)-1)/opts.scale)*0.03; % extremity at 3% of the bottom edge of the frame
            plot([scaleBarXmin; scaleBarXmax], [scaleBarY; scaleBarY], '-k', 'LineWidth', 3) % change the color of the bar as needed. '-k' for black, '-w' for white
            text(scaleBarXmin,(size(I,1)-1)/opts.scale*0.05,sprintf('%i cm',scalebarLength*100),'Color',[0 0 0]) % change the text and its color as needed

            %---Reference vector---
            p1 = [((size(I,2)-1)/opts.scale)*0.03 ((size(I,1)-1)/opts.scale)*0.03]; % First Point at 3% of max y axis and x axis (x and y of origin of scale vector)
            quiver(p1(1),p1(2),0,opts.forceScale*arrowLength,0,'MaxHeadSize',0.7,'LineWidth', 2,'Color', 'k') % reference vector scale
            text(((size(I,2)-1)/opts.scale)*0.035,((size(I,1)-1)/opts.scale)*0.04,sprintf('%.3f mN cm^{-1}',arrowLength),'Color',[0 0 0])
        else
            %---Scalebar---
            % if the axes are turned on it may not be necessary. But if ensures the
            % figure is properly formatted all the time
            % scalebarLength = 0.02; % scale (2cm), change as needed
            scaleBarXmax = ((size(I,2)-1)/opts.scale)*0.97; % extremity at 97% of the right edge of the frame
            scaleBarXmin = scaleBarXmax-scalebarLength; %  
            scaleBarY = ((size(I,1)-1)/opts.scale)*0.03; % extremity at 3% of the bottom edge of the frame
            plot([scaleBarXmin; scaleBarXmax], [scaleBarY; scaleBarY], '-w', 'LineWidth', 3) % change the color of the bar as needed. '-k' for black, '-w' for white
            text(scaleBarXmin,(size(I,1)-1)/opts.scale*0.05,sprintf('%i cm',scalebarLength*100),'Color',[1 1 1]) % change the text and its color as needed

            %---Reference vector---
            p1 = [((size(I,2)-1)/opts.scale)*0.03 ((size(I,1)-1)/opts.scale)*0.03]; % First Point at 3% of max y axis and x axis (x and y of origin of scale vector)
            quiver(p1(1),p1(2),0,opts.forceScale*arrowLength,0,'MaxHeadSize',0.7,'LineWidth', 2,'Color', 'w') % reference vector scale
            text(((size(I,2)-1)/opts.scale)*0.035,((size(I,1)-1)/opts.scale)*0.04,sprintf('%.3f mN cm^{-1}',arrowLength),'Color',[1 1 1])
       end

    title(sprintf('Axial force - frame %i',fr))
    axis([0,(size(I,2)-1)/opts.scale, 0,(size(I,1)-1)/opts.scale])
    axis off
    axis equal
    formatFigure

subplot(1,3,2)
    hold on
            if opts.useRawImage == 0 
               patch(boundx,boundy,'k')
            else
               imagesc([0,(size(I,2)-1)/opts.scale],[0,(size(I,1)-1)/opts.scale],flipud(I)); colormap('gray') % plot image 
            end

    quiver(boundx2(indleftpull,:),boundy2(indleftpull,:),surfpress(indleftpull,:).*Lx(indleftpull,:).*darea(indleftpull,:)*opts.forceScale,surfpress(indleftpull,:).*Ly(indleftpull,:).*darea(indleftpull,:)*opts.forceScale,0,'Color','#562689','LineWidth',arrowThickness,'MaxHeadSize',arrowHead);   % original color was [0 1 1]
    quiver(boundx2(indleftpush,:),boundy2(indleftpush,:),surfpress(indleftpush,:).*Lx(indleftpush,:).*darea(indleftpush,:)*opts.forceScale,surfpress(indleftpush,:).*Ly(indleftpush,:).*darea(indleftpush,:)*opts.forceScale,0,'Color','#b35807','LineWidth',arrowThickness,'MaxHeadSize',arrowHead);
    quiver(boundx2(indrightpull,:),boundy2(indrightpull,:),surfpress(indrightpull,:).*Lx(indrightpull,:).*darea(indrightpull,:)*opts.forceScale,surfpress(indrightpull,:).*Ly(indrightpull,:).*darea(indrightpull,:)*opts.forceScale,0,'Color','#562689','LineWidth',arrowThickness,'MaxHeadSize',arrowHead);
    quiver(boundx2(indrightpush,:),boundy2(indrightpush,:),surfpress(indrightpush,:).*Lx(indrightpush,:).*darea(indrightpush,:)*opts.forceScale,surfpress(indrightpush,:).*Ly(indrightpush,:).*darea(indrightpush,:)*opts.forceScale,0,'Color','#b35807','LineWidth',arrowThickness,'MaxHeadSize',arrowHead);

        if opts.useRawImage == 0 
            %---Scalebar---
            % if the axes are turned on it may not be necessary. But if ensures the
            % figure is properly formatted all the time
            % scalebarLength = 0.02; % scale (2cm), change as needed
            scaleBarXmax = ((size(I,2)-1)/opts.scale)*0.97; % extremity at 97% of the right edge of the frame
            scaleBarXmin = scaleBarXmax-scalebarLength; %  
            scaleBarY = ((size(I,1)-1)/opts.scale)*0.03; % extremity at 3% of the bottom edge of the frame
            plot([scaleBarXmin; scaleBarXmax], [scaleBarY; scaleBarY], '-k', 'LineWidth', 3) % change the color of the bar as needed. '-k' for black, '-w' for white
            text(scaleBarXmin,(size(I,1)-1)/opts.scale*0.05,sprintf('%i cm',scalebarLength*100),'Color',[0 0 0]) % change the text and its color as needed

            %---Reference vector---
            p1 = [((size(I,2)-1)/opts.scale)*0.03 ((size(I,1)-1)/opts.scale)*0.03]; % First Point at 3% of max y axis and x axis (x and y of origin of scale vector)
            quiver(p1(1),p1(2),0,opts.forceScale*arrowLength,0,'MaxHeadSize',0.7,'LineWidth', 2,'Color', 'k') % reference vector scale
            text(((size(I,2)-1)/opts.scale)*0.035,((size(I,1)-1)/opts.scale)*0.04,sprintf('%.3f mN cm^{-1}',arrowLength),'Color',[0 0 0])
        else
            %---Scalebar---
            % if the axes are turned on it may not be necessary. But if ensures the
            % figure is properly formatted all the time
            % scalebarLength = 0.02; % scale (2cm), change as needed
            scaleBarXmax = ((size(I,2)-1)/opts.scale)*0.97; % extremity at 97% of the right edge of the frame
            scaleBarXmin = scaleBarXmax-scalebarLength; %  
            scaleBarY = ((size(I,1)-1)/opts.scale)*0.03; % extremity at 3% of the bottom edge of the frame
            plot([scaleBarXmin; scaleBarXmax], [scaleBarY; scaleBarY], '-w', 'LineWidth', 3) % change the color of the bar as needed. '-k' for black, '-w' for white
            text(scaleBarXmin,(size(I,1)-1)/opts.scale*0.05,sprintf('%i cm',scalebarLength*100),'Color',[1 1 1]) % change the text and its color as needed

            %---Reference vector---
            p1 = [((size(I,2)-1)/opts.scale)*0.03 ((size(I,1)-1)/opts.scale)*0.03]; % First Point at 3% of max y axis and x axis (x and y of origin of scale vector)
            quiver(p1(1),p1(2),0,opts.forceScale*arrowLength,0,'MaxHeadSize',0.7,'LineWidth', 2,'Color', 'w') % reference vector scale
            text(((size(I,2)-1)/opts.scale)*0.035,((size(I,1)-1)/opts.scale)*0.04,sprintf('%.3f mN cm^{-1}',arrowLength),'Color',[1 1 1])
        end

    title(sprintf('Lateral force - frame %i',fr))
    axis([0,(size(I,2)-1)/opts.scale, 0,(size(I,1)-1)/opts.scale])
    axis off
    axis equal
    formatFigure

subplot(1,3,3)
    hold on
        if opts.useRawImage == 0 
            patch(boundx,boundy,'k')
        else
            imagesc([0,(size(I,2)-1)/opts.scale],[0,(size(I,1)-1)/opts.scale],flipud(I)); colormap('gray') % plot image

        end

    quiver(boundx2(indpospull,:),boundy2(indpospull,:),opts.forceScale*surfpress(indpospull,:).*surfunitnormx(indpospull,:).*darea(indpospull,:),opts.forceScale*surfpress(indpospull,:).*surfunitnormy(indpospull,:).*darea(indpospull,:),0,'Color',[0.45 0.68 0.8],'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);
    quiver(boundx2(indpospush,:),boundy2(indpospush,:),opts.forceScale*surfpress(indpospush,:).*surfunitnormx(indpospush,:).*darea(indpospush,:),opts.forceScale*surfpress(indpospush,:).*surfunitnormy(indpospush,:).*darea(indpospush,:),0,'Color',[0.8 0.5 0.5],'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);
    quiver(boundx2(indnegpull,:),boundy2(indnegpull,:),opts.forceScale*surfpress(indnegpull,:).*surfunitnormx(indnegpull,:).*darea(indnegpull,:),opts.forceScale*surfpress(indnegpull,:).*surfunitnormy(indnegpull,:).*darea(indnegpull,:),0,'Color',[0 0.2 0.5],'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);
    quiver(boundx2(indnegpush,:),boundy2(indnegpush,:),opts.forceScale*surfpress(indnegpush,:).*surfunitnormx(indnegpush,:).*darea(indnegpush,:),opts.forceScale*surfpress(indnegpush,:).*surfunitnormy(indnegpush,:).*darea(indnegpush,:),0,'Color',[0.5 0 0.1],'LineWidth',arrowThickness,'MaxHeadSize',arrowHead)

        if opts.useRawImage == 0 
            %---Scalebar---
            % if the axes are turned on it may not be necessary. But if ensures the
            % figure is properly formatted all the time
            % scalebarLength = 0.02; % scale (2cm), change as needed
            scaleBarXmax = ((size(I,2)-1)/opts.scale)*0.97; % extremity at 97% of the right edge of the frame
            scaleBarXmin = scaleBarXmax-scalebarLength; %  
            scaleBarY = ((size(I,1)-1)/opts.scale)*0.03; % extremity at 3% of the bottom edge of the frame
            plot([scaleBarXmin; scaleBarXmax], [scaleBarY; scaleBarY], '-k', 'LineWidth', 3) % change the color of the bar as needed. '-k' for black, '-w' for white
            text(scaleBarXmin,(size(I,1)-1)/opts.scale*0.05,sprintf('%i cm',scalebarLength*100),'Color',[0 0 0]) % change the text and its color as needed

            %---Reference vector---
            p1 = [((size(I,2)-1)/opts.scale)*0.03 ((size(I,1)-1)/opts.scale)*0.03]; % First Point at 3% of max y axis and x axis (x and y of origin of scale vector)
            quiver(p1(1),p1(2),0,opts.forceScale*arrowLength,0,'MaxHeadSize',0.7,'LineWidth', 2,'Color', 'k') % reference vector scale
            text(((size(I,2)-1)/opts.scale)*0.035,((size(I,1)-1)/opts.scale)*0.04,sprintf('%.3f mN cm^{-1}',arrowLength),'Color',[0 0 0])
        else
            %---Scalebar---
            % if the axes are turned on it may not be necessary. But if ensures the
            % figure is properly formatted all the time
            % scalebarLength = 0.02; % scale (2cm), change as needed
            scaleBarXmax = ((size(I,2)-1)/opts.scale)*0.97; % extremity at 97% of the right edge of the frame
            scaleBarXmin = scaleBarXmax-scalebarLength; %  
            scaleBarY = ((size(I,1)-1)/opts.scale)*0.03; % extremity at 3% of the bottom edge of the frame
            plot([scaleBarXmin; scaleBarXmax], [scaleBarY; scaleBarY], '-w', 'LineWidth', 3) % change the color of the bar as needed. '-k' for black, '-w' for white
            text(scaleBarXmin,(size(I,1)-1)/opts.scale*0.05,sprintf('%i cm',scalebarLength*100),'Color',[1 1 1]) % change the text and its color as needed

            %---Reference vector---
            p1 = [((size(I,2)-1)/opts.scale)*0.03 ((size(I,1)-1)/opts.scale)*0.03]; % First Point at 3% of max y axis and x axis (x and y of origin of scale vector)
            quiver(p1(1),p1(2),0,opts.forceScale*arrowLength,0,'MaxHeadSize',0.7,'LineWidth', 2,'Color', 'w') % reference vector scale
            text(((size(I,2)-1)/opts.scale)*0.035,((size(I,1)-1)/opts.scale)*0.04,sprintf('%.3f mN cm^{-1}',arrowLength),'Color',[1 1 1])
        end

        title(sprintf('Force vectors - frame %i',fr))
        axis([0,(size(I,2)-1)/opts.scale, 0,(size(I,1)-1)/opts.scale])
        axis off
        axis equal
        formatFigure
end

if firstLast(1)~=firstLast(2)
    if opts.exportFigures == 1
        % In general, it is fine just exporting simple jpeg figures. However, I
        % always export .svg files that are fully editable in Illustrator
        % because they are vector files (lossless format). Select the option
        % that fits your needs
        % Three export options: 1) Export only jpg; 2) export only .svg; 3) export both .jpg and .svg.
        if opts.exportPlotFormat == 1 % Export as .jpg
            FilepathJPG = fullfile(PathForcesFigures,sprintf('force%04d',fr));
            print('-djpeg', FilepathJPG)
    
        elseif opts.exportPlotFormat == 2 % Exports as .svg (preferred)
            FilepathSVG = fullfile(PathForcesFigures,sprintf('force%04d',fr));
            print('-dsvg', '-vector',FilepathSVG) 
    
        elseif opts.exportPlotFormat ==3 % exporting .jpg and .svg together will take longer
            FilepathJPG = fullfile(PathForcesFigures,sprintf('force%04d',fr));
            print('-djpeg', FilepathJPG)
            FilepathSVG = fullfile(PathForcesFigures,sprintf('force%04d',fr));
            print('-dsvg', '-vector',FilepathSVG) 
    
        else % in case no number is provided
            FilepathJPG = fullfile(PathForcesFigures,sprintf('force%04d',fr));
            print('-djpeg', FilepathJPG)
        end
    end
end

if firstLast(1)~=firstLast(2)
   pause(0.1) % leaves time to visualize the plot. WARNING: increases the duration of the export
   hold off
   clf % clear figure to conserve memory
end
%% Calculate force acting on the body
% Axial force
totypospull = sum(darea(indpospull,:).*surfpress(indpospull,:).*surfunitnormy2(indpospull,:));  % total forward pull force in this frame
totypospush = sum(darea(indpospush,:).*surfpress(indpospush,:).*surfunitnormy2(indpospush,:));  % total forward push force in this frame
totynegpull = sum(darea(indnegpull,:).*surfpress(indnegpull,:).*surfunitnormy2(indnegpull,:));  % total backward pull force in this frame
totynegpush = sum(darea(indnegpush,:).*surfpress(indnegpush,:).*surfunitnormy2(indnegpush,:));  % total backward push force in this frame
netyforce = nansum(darea(:,:).*surfpress(:,:).*surfunitnormy2(:,:)); % net force in this frame

% Lateral force
totxleftpull = sum(darea(indleftpull,:).*surfpress(indleftpull,:).*surfunitnormx2(indleftpull,:));  % total left pull force in this frame
totxleftpush = sum(darea(indleftpush,:).*surfpress(indleftpush,:).*surfunitnormx2(indleftpush,:));  % total left push force in this frame
totxrightpull = sum(darea(indrightpull,:).*surfpress(indrightpull,:).*surfunitnormx2(indrightpull,:));  % total right pull force in this frame
totxrightpush = sum(darea(indrightpush,:).*surfpress(indrightpush,:).*surfunitnormx2(indrightpush,:));  % total right push force in this frame
netxforce = nansum(darea(:,:).*surfpress(:,:).*surfunitnormx2(:,:)); % net force in this frame       

%% Compute power and efficiency
if fr > 1 % cannot compute power of first frame because we need to know the position of the body in frame n-1
    % Import previous blanking outline file
    prevblank = importdata(quickfilepath(D.outlines(fr-opts.increment))); %D.outlines(fr-opts.increment); % importdata(quickfilepath(D_outlines(fr-opts.increment))); % import previous outline blanking coordinates
    prevblank = curvspace(prevblank,opts.npoints);
    prevblank = alignPoints(prevblank); % always align array with anteriormost point of animal as first element 
    prevboundx = prevblank(:,1)/1000;  % convert x coordinate from mm to meters
    prevboundy = prevblank(:,2)/1000;  % convert y coordinate from mm to meters
   
    %calculate local speed of body surface
    surfvelu = abs((boundx-prevboundx)./deltat); % x-direction
    surfvelv = abs((boundy-prevboundy)./deltat); % y-direction
    
    % calculate local force of animal on fluid
    localforcex = darea.*surfpress.*surfunitnormx2; % x-direction
    localforcey = darea.*surfpress.*surfunitnormy2; % y-direction
    
    % Calculate local power exerted by animal
    locallatpower = abs(surfvelu.*localforcex);   % in lateral direction (assuming animal swims in Y-DIRECTION)
    localaxipower = abs(surfvelv.*localforcey);   % in axial direction (assuming animal swims in Y-DIRECTION)
    
    % Calculate total power
    latpower = nansum(locallatpower,1); % lateral 
    axipower = nansum(localaxipower,1); % axial
    
    % Calculate total power exerted due to negative/positive suction/push pressure
    latpullpower = nansum(locallatpower(find(surfpress<0),1)); % total lateral power exerted due to low pressure suction
    latpushpower = nansum(locallatpower(find(surfpress>0),1)); % total lateral power exerted due to high pressure pushing
    axipullpower = nansum(localaxipower(find(surfpress<0),1)); % total axial power exerted due to low pressure suction
    axipushpower = nansum(localaxipower(find(surfpress>0),1)); % total axial power exerted due to high pressure pushing
    
    % Calculate fish swiming efficiency
    %swimmingSpeed = 0.060675; % average animal swimming speed in meter
    FroudeEfficiency = ((totypospull+totypospush)*opts.swimmingSpeed)/(((totypospull+totypospush)*opts.swimmingSpeed)+latpower);
    
else
    latpower = 0;
    axipower = 0;
    localforcex = darea.*surfpress.*surfunitnormx2; % x-direction
    localforcey = darea.*surfpress.*surfunitnormy2; % y-direction
    locallatpower = zeros(opts.npoints,1);
    localaxipower = zeros(opts.npoints,1);
    latpullpower = 0;
    latpushpower = 0;
    axipullpower = 0;
    axipushpower = 0;
    FroudeEfficiency = 0;
end

% Store summary force and power data for each frame in a new variable
% frame #, total pull thrust, total push thrust, total pull drag, total
% push drag, net force, total axial power, total lateral power, total axial
% pull power, total axial push power, total lateral pull power,  total lat
% push power, Froude efficiency
summaryData(advanceFr,:) = [fr,totypospull,totypospush,totynegpull,totynegpush,netyforce,axipower,latpower,axipullpower,axipushpower,latpullpower,latpushpower,FroudeEfficiency];

% Store local force and power along the body in cell array
% x coordinates (m), y coordinates (m), local lateral force component (N m-1), local axial
% force component (N m-1), local lateral power (W m-1), local axial power (W m-1)
localForcePowerData = [boundx2,boundy2,localforcex,localforcey,locallatpower,localaxipower]; % concatenate data
allLocalForceData(advanceFr) = {localForcePowerData};

% Export local force and power along the body (1 file per frame)
if firstLast(1)~=firstLast(2)
filename = sprintf('local_%04d.csv',fr); % file name and .csv file type
% dlmwrite(fullfile(PathLocalForces,filename),localForcePowerData);
writematrix(localForcePowerData,fullfile(PathLocalForces,filename))
end

advanceFr = advanceFr+1; % advance index of summaryData by 1    
end

if firstLast(1)~=firstLast(2)
% Export final summary table
dataTable = array2table(summaryData);
dataTable.Properties.VariableNames = {'frame_#','total_pull_thrust_N.m-1','total_push_thrust_N.m-1','total_pull_drag_N.m-1','total_push_drag_N.m-1','net_force_N.m-1','total_axial_power_W.m-1','total_lateral_power_W.m-1','total_axial_pull_power_W.m-1','total_axial_push_power_W.m-1','total_lateral_pull_power_W.m-1','total_lateral_push_power_W.m-1','Froude efficiency'}; % set the name of each variable
% dataTable_out = fullfile(PathForces,'Forces and power summary.csv');
writetable(dataTable,fullfile(path.forces,'Forces and power summary.csv')) % write table in the file format defined in the line above (.csv)
fprintf('Summary table exported\n');
fprintf('Local force and power data exported\n');
close % close empty remaining figure
end

end




