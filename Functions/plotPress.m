function plotPress(pressData,firstLast,opts,path,D) % pressData contains the meshgrid X,Y,P for all the frames
if firstLast(1) ~= firstLast(2) % if the first and last frame are different, then we are in the-all plot export mode. We thus need to generate the export paths for the figures.
    % In general, it is fine just exporting simple jpeg figures. However, I
    % always export .svg files that are fully editable in Illustrator
    % because they are vector files (lossless format). Select the option
    % that fits your needs
    % Three export options: 1) Export only jpg; 2) export only .svg; 3) export both .jpg and .svg.
    if opts.exportPlot == 1
     filepath_jpg = fullfile(path.main,'Figures','Pressure','pressMask-jpg');
     mkdir(filepath_jpg)
    elseif opts.exportPlot == 2
     filepath_svg = fullfile(path.main,'Figures','Pressure','pressMask-svg');
     mkdir(filepath_svg)
    elseif opts.exportPlot ==3
     filepath_jpg = fullfile(path.main,'Figures','Pressure','pressMask-jpg');
     filepath_svg = fullfile(path.main,'Figures','Pressure','pressMask-svg');
     mkdir(filepath_svg)
     mkdir(filepath_jpg)
    else % in case no number is provided
     filepath_jpg = fullfile(path.main,'Figures','Pressure','pressMask-jpg');
     mkdir(filepath_jpg)
    end
end

% Set custom pressure map
press_colormap = customcolormap_preset('orange-white-purple'); % define the pressure colormap (use custom colors to enhance the visual clarity of the data)
hLines = 201; % number of isolines; change according to how smooth you want the colorbar and the plot to look like. Over 200 color steps will slow down plotting.

% Increase the resolution of the pressure fields to ensure a smooth cutout of the field by the mask.
% Option to change the resolution in two different ways. 1 = set the increment between vectors by defining the actual physical spacing in meter; 2) Increase the resolution by a factor of n. The higher the resolution, the longer it will tekae to interpolate the fields.

if opts.Resolution == 1
    % opts.smoothness sets the spacing between vectors on your scale (in meters).
    [Xq,Yq] = meshgrid(min(pressData.X(:)):opts.smoothness:max(pressData.X(:)),min(pressData.Y(:)):opts.smoothness:max(pressData.Y(:))); % makes Y and Y data grid with higher resolution

elseif opts.Resolution == 2
    % opts.smoothness increase the resolution by a factor of n
    distVec = (pressData.X(1,end)-pressData.X(1,end-1))/opts.smoothness; % calculate the distance between vectors for the correspinding smoothness factor
    [Xq,Yq] = meshgrid(min(pressData.X(:)):distVec:max(pressData.X(:)),min(pressData.Y(:)):distVec:max(pressData.Y(:))); % makes Y and Y data grid with higher resolution

else
    Xq = veloData.X;
    Yq = veloData.Y;
end

%-----Plot figure-----
fig = figure;
set(gcf, 'WindowState', 'maximized');

% plotting loop (during testing, because the first and last frames are the
% same, the loop will only process one frame)
if firstLast(1) ~= firstLast(2) % if the first and last frame are different, then we are in the-all plot export mode.
fprintf('Plotting pressure fields...');
end

for ff=firstLast(1):firstLast(2) % run loop from first frame to last frame listed in the second function input)
    hold on % to overlay data in the figure. keep it in the loop because hold off is added at the end to avoid potential issues if checking multiple random frames during testing.

%---load presure field for the current frame---
    press = pressData.P(:,:,ff);
    Pq = interp2(pressData.X,pressData.Y,press,Xq,Yq); % interpolates pressure data from original XY grid to higher resolution grid

%---use outlines---
    if opts.outline == 1 % use actual xy outlines
    outlines = importdata(quickfilepath(D.outlines(ff))); % outlines should have been scaled to meters
    
    elseif opts.outline == 2 % use BW masks
    % Find the outlines from the BW masks
    BWtemp = flipud(importdata(quickfilepath(D.BW(ff)))); % import mask
    BWdiag = bwmorph(BWtemp,'diag'); % diagonal fill
    BWtempDilate = bwmorph(BWdiag,'thicken',1); % dilate mask by 1px so the outline falls exactly on the edge of the mask
    if sum(BWtempDilate(:))>0
    B = bwboundaries(BWtempDilate,8);
    outlinesRaw = cell(length(B)*2-1,1);
    a = 1;
    for ii = 1:length(B)
    outlinesRaw{a,1} = B{ii};
    outlinesRaw{a+1,1} = [NaN,NaN];
    a = a+2;
    end
    outlines = cell2mat(outlinesRaw)/opts.scale; % final matrix containing the outlines with holes;
    
    else % if no masking is required (not recommended)
    outlines = [];
    end
    end

%--- Delete pressure within the mask---
    if opts.outline == 1
        if isempty(outlines) == 0
            [in] = inpolygon(Xq,Yq,outlines(:,1),outlines(:,2)); % finds XY data in high-res grid that falls in-on mask
            in = (in-1).^2; % transforms 1 to 0 for conversion to NaN, and 0 to 1
            in(in == 0) = NaN;
            Pq = Pq.*in;
        end
    elseif opts.outline == 2
        if isempty(outlines) == 0
            [in] = inpolygon(Xq,Yq,outlines(:,2),outlines(:,1)); % finds XY data in high-res grid that falls in-on mask
            in = (in-1).^2; % transforms 1 to 0 for conversion to NaN, and 0 to 1
            in(in == 0) = NaN;
            Pq = Pq.*in;
        end
    end

%---Import the raw image---
I = importdata(quickfilepath(D.images(ff))); % import corresponding image

%---Overlay the pressure fields---
opacity = 1; % opacity factor to control the ransparency of the pressure fields. 0 = completely transparent, 1 = completely opaque. Useful when asked to show the particle field.
imagesc([0,(size(I,2)-1)/opts.scale],[0,(size(I,1)-1)/opts.scale],flipud(I)); colormap('gray'); % offset outline by the equivalent of 1px because bwboundaries starts the image at (x,y)=1 and outputs the outline at (x,y)=1.5. The outline should be set at the center of the poixel, thus (x,y)=0.5 when the origin of the BW mask is set to 0. Set the scale for the outline, in m. 
contourf(Xq,Yq,Pq,hLines,'edgecolor','none','FaceAlpha',opacity);   % Plot pressure fields; change the opacity if needed (from 0 to 1)

%---prepare polyshape for outline--- 
    if opts.outline == 1
        if isempty(outlines) == 0
    pgon = polyshape(outlines(:,1),outlines(:,2));
        end
    elseif opts.outline == 2
        if isempty(outlines) == 0
    pgon = polyshape(outlines(:,2),outlines(:,1));
        end
    end

%---plot actual mask if needed---
    if opts.plotMask == 1
        if isempty(outlines) == 0
    plot(pgon,'FaceColor','white','FaceAlpha',1,'EdgeColor','black') % default color is white infill and black edge. The Opacity of the mask is controlled by 'FaceAlpha'
        end
    end

%---Set the figure axes length---
axis equal
if opts.tightPlot == 1
   % Change parameters as needed
    axisXmin = 0.004; %was 0.040 (in m)
    axisYmin = 0.04; %was 0.005 (in m)
    axisXmax = 0.126; % was 0.05
    axisYmax = 0.12; % was 0.05
    axis([axisXmin axisXmax axisYmin axisYmax])

else % plot the whole frame
    axisXmin = 0; 
    axisYmin = 0; 
    axisXmax = max(unique(pressData.X));
    axisYmax = max(unique(pressData.Y));

    % Plot properties - not necessary when axes are turned off
%     xlabel('distance (m)')
%     ylabel('distance (m)')
end

%---Scalebar---
% if the axes are turned on it may not be necessary. But if ensures the
% figure is properly formatted all the time

scalebarLength = 0.02; % scale (2cm), change as needed
scaleBarXmax = axisXmax*0.97; % extremity at 97% of the right edge of the frame
scaleBarXmin = scaleBarXmax-scalebarLength; %  
scaleBarY = axisYmax*0.03; % extremity at 5% of the bottom edge of the frame
plot([scaleBarXmin; scaleBarXmax], [scaleBarY; scaleBarY], '-k', 'LineWidth', 3) % change the color of the bar as needed. '-k' for black, '-w' for white
text(scaleBarXmin,axisYmax*0.05,sprintf('%i cm',scalebarLength*100),'Color',[0 0 0]) % change the text and its color as needed


%---Colormap properties---
colormap(press_colormap)
hcb = colorbar;
clim([-opts.pressMag opts.pressMag]) % set the min and max vorticity magnitude in s-1
hcb.Label.String = 'pressure (Pa)';
hcb.Label.FontSize = 12;
hcb.Label.FontName = 'Arial';

%---set the figure axes---
% keep at the end to ensure that vectors extending outside the domain don't
% deform the figure
xlim([axisXmin axisXmax])
ylim([axisYmin axisYmax])
axis off % turn the axis off for clarity

formatFigure
hold off

if firstLast(1) ~= firstLast(2) % if the first and last frame are different, then we are in the all-plot export mode. We thus need to generate the export paths for the figures.
    % In general, it is fine just exporting simple jpeg figures. However, I
    % always export .svg files that are fully editable in Illustrator
    % because they are vector files (lossless format). Select the option
    % that fits your needs
    % Three export options: 1) Export only jpg; 2) export only .svg; 3) export both .jpg and .svg.
    if opts.exportPlot == 1 % Export as .jpg
        FilepathJPG = fullfile(filepath_jpg,sprintf('press%04d',ff));
        print('-djpeg', FilepathJPG)

    elseif opts.exportPlot == 2 % Exports as .svg (preferred)
        FilepathSVG = fullfile(filepath_svg,sprintf('press%04d',ff));
        print('-dsvg', '-vector',FilepathSVG) 

    elseif opts.exportPlot ==3 % exporting .jpg and .svg together will take longer
        FilepathJPG = fullfile(filepath_jpg,sprintf('press%04d',ff));
        print('-djpeg', FilepathJPG)
        FilepathSVG = fullfile(filepath_svg,sprintf('press%04d',ff));
        print('-dsvg', '-vector',FilepathSVG) 

    else % in case no number is provided
        FilepathJPG = fullfile(filepath_jpg,sprintf('press%04d',ff));
        print('-djpeg', FilepathJPG)
    end
% disp(ff) 
clf % clear figure so we don't overlap sevral together

if firstLast(1) ~= firstLast(2) % if the first and last frame are different, then we are in the-all plot export mode.
progressCount(ff,firstLast(2)-firstLast(1)+1); % display export progress

if ff == firstLast(2) % when the loop reaches the last frame during export, simply close the figure and confirm the export was successful
    close(fig)
    fprintf('Done\n');
    beep on; beep

end
end

end

end