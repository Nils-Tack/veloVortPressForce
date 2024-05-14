function plotFrameWithOutline(I,outline,opts)
    figure; hold on
    imagesc([0,(size(I,2)-1)/opts.scale],[0,(size(I,1)-1)/opts.scale],flipud(I)); colormap('gray') % plot image
    
    if size(outline,1)>1 % only plot if the outline file contains at least 2 points; if the outline file is empty because no mask was created or contyains only a default 0,0 point, then simply do not plot
    plot(outline(:,1), outline(:,2), 'r','LineWidth',2);
    else
        axisXmax = ((size(I,2)-1)/opts.scale)*0.03; % origin of the text at 3% of the left edge of the frame
        text(axisXmax,axisXmax,'No outline data for this image, move to the next frame','Color',[1 0 0]) % change the text and its color as needed
    end

    % Label the axes to remember the scale
    xlabel('distance (m)')
    ylabel('distance (m)')
    axis equal

end