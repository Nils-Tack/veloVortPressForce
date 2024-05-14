function outlineFixed = fixOutlinePoints(outline,I,zoomAxisLim,opts)

% correct brightness of the image to incease visibility
I = I+70;

% Plot the figure with the outline
figure; hold on
imagesc([0,(size(I,2)-1)/opts.scale],[0,(size(I,1)-1)/opts.scale],flipud(I)); colormap('gray') % plot image
plot(outline(:,1), outline(:,2), '.r','MarkerSize',10);
axis('equal')
axis(zoomAxisLim)
figSize(1.2,1.2,gcf)

title("Repair desired section of outline. Hit 'return' to finish.")
[xd,yd] = ginputc('ShowPoints',true,'Color','r'); % Use red crosshairs; select as many point as needed and hit return when done

% Find nearest neighbors
nn = knnsearch(outline,[xd yd]);

% Get new section
newsect = [xd,yd];

% Replace ends with nns
newsect(1,:) = outline(nn(1),:); % to perform a perfect blend of the outline with the new section (first point)
newsect(end,:) = outline(nn(end),:); % to perform a perfect blend of the outline with the new section (last point)

% Interpolate
newsect = interp2path(unique(newsect,'stable','row'),'pchip',1);

% Replace points
if nn(1) < nn(end)
% Split mask into two sections, excluding new section
mask_start = outline(1:nn(1),:);
mask_end = outline(nn(end):end,:);

% Fill in section
outlineFixed = [mask_start; newsect(2:end-1,:); mask_end]; %original 
    
else
mask_start = outline(1:nn(end),:);% was 1
mask_end = outline(nn(1):end,:);% was end

% Fill in section
outlineFixed = [mask_start; flipud(newsect(2:end-1,:)); mask_end]; %original 

end

outlineFixed = unique(outlineFixed,'stable','row');

end