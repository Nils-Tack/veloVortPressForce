function outline = makeOutlines(BW,opts)

BW = flipud(BW); % flip the image up/down to maintain the origin of the y axis in the bottom left corner
BWdiag = bwmorph(BW,'diag'); % diagonal fill to smooth the outline and avoid crenated edges caused by the pixels
BWtempDilate = bwmorph(BWdiag,'thicken',1); % dilate mask by 1px so the coordinate of the outline point fall exactly on the edge of the mask

if sum(BWtempDilate(:))>0 % necessary to check that there is a mask in the image
BWedge = bwboundaries(BWtempDilate,'noholes'); % Find edge of the subject

    if opts.multiOutlines
        outlinesRaw = cell(length(BWedge)*2-1,1);
        a = 1;
            for ii = 1:length(BWedge)
            outlinesRaw{a,1} = BWedge{ii};
            outlinesRaw{a+1,1} = [NaN,NaN];
            a = a+2;
            end
        outline = fliplr((cell2mat(outlinesRaw)/opts.scale)-(1/opts.scale)); % final matrix containing the edge of multiple outlines
    else
        outline = [(BWedge{1}(:,2)/opts.scale)-(1/opts.scale),(BWedge{1}(:,1)/opts.scale)-(1/opts.scale)]; % offset outline by the equivalent of 1px because bwboundaries starts the image at (x,y)=1 and outputs the outline at (x,y)=1.5. The outline should be set at the center of the pixel, thus (x,y)=0.5 when the origin of the BW mask is set to 0. Set the scale for the outline, in m. 
        N = 500; % number of interpolated points
        outline = curvspace([outline(:,1) outline(:,2)],N);

    end
else
    outline = [];
end

% Plot the mask and outline
hold on
imagesc([0,(size(BW,2)-1)/opts.scale],[0,(size(BW,1)-1)/opts.scale],BWtempDilate); colormap('gray') % plot mask using BW image directly
if sum(BW(:))>0 % check that there is a mask in the figure, and thus a mask.
plot(outline(:,1), outline(:,2), 'r','LineWidth',2);
end

% Label the axes to remember the scale
xlabel('distance (m)')
ylabel('distance (m)')
axis([0 (size(BW,2)-1)/opts.scale 0 (size(BW,1)-1)/opts.scale])
axis equal


end