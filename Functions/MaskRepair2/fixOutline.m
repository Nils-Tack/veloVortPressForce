function outlineFixed = fixOutline(I,outline, opts)
% Select repair location and start repair on outline
% Options
zoomWindowSize = 0.05; % Size of the zoomed field of view on the image in meter (square window of size zoomWindowSize x zoomWindowSize)

% Show outline
figure; hold on % Create new figure with the outline that will be modified
imagesc([0,(size(I,2)-1)/opts.scale],[0,(size(I,1)-1)/opts.scale],flipud(I)); colormap('gray') % plot image
plot(outline(:,1), outline(:,2), '.r','MarkerSize',10);
axis('equal')

figSize(1.2,1.2,gcf)
title('Click near section to repair.')

% Get axes around section to repair
[xax, yax] = ginput(1); % select a point around where to fix the outline
zoomAxisLim = [xax-zoomWindowSize/2, xax+zoomWindowSize/2, yax-zoomWindowSize/2, yax+zoomWindowSize/2]; % find the limit of the axes to display a zoomed section of the image and outline

% Replace points
fprintf('Replacing points.')
outlineFixed = fixOutlinePoints(outline,I,zoomAxisLim,opts);

end