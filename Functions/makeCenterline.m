function centerline = makeCenterline(BW,opts)

%---Compute the outline from the BW---
BWdiag = bwmorph(BW,'diag'); % diagonal fill to smooth the outline and avoid crenated edges caused by the pixels
BWtempDilate = bwmorph(BWdiag,'thicken',1); % dilate mask by 1px so the coordinate of the outline point fall exactly on the edge of the mask
BWedge = bwboundaries(BWtempDilate,'noholes'); % Find edge of the subject
outline = [BWedge{1}(:,2),BWedge{1}(:,1)]; % offset outline by the equivalent of 1px because bwboundaries starts the image at (x,y)=1 and outputs the outline at (x,y)=1.5. The outline should be set at the center of the pixel, thus (x,y)=0.5 when the origin of the BW mask is set to 0. Set the scale for the outline, in m. 


%---Skeletonize the mask---
BWskel = bwskel(BW,'MinBranchLength',180); % change value to eliminate branching if needed. The higher the value (in px) the longer the branches that are removed
BWskel = bwmorph(BWskel,'bridge'); % Just in case the centerline gets interupted by an empty pixel at the branching upon discarding branch, add connectivity.

%---Find coordinates of binary centerline points---
[row,col] = find(BWskel);

%---Order the indices of the points along centerline in sequence---
    % Find the ends of the centerline
    BWends = bwmorph(BWskel,'endpoints');
    [r,c] = find(BWends); % find both ends but select only one
    
    % Find the index of the initial point in the row-col matrix using the
    % extremities identified from the skeleton
    tempIndex = find(row == r(2) & col == c(2));
    
    % circshift the row and col to have the sequence start exactly at the first point
    row = circshift(row,-tempIndex+1);
    col = circshift(col,-tempIndex+1);

    % Uses closest-neighbor method to organize the indices of each point in sequence theough propagation from the extremity identified above. 
    % This method is appropriate in most cases, at least in the context of clean fish centerlines
    % Set up loop variables. The first coordinate will be the first coordinate in the new, sorted order.
%     yIdx = [1,nan(size(row)-1)]; 
    xyTemp = [col(:), row(:)]; 
    xyTemp(1,:) = NaN; 
    idx = [1, nan(1, numel(col)-1)]; 
    counter = 0; 
    % Loop through each coordinate and find its nearest neightbor,
    % Then eliminate that coordinate from being chosen again, and
    % then start the loop again finding thenearest neighbor to 
    % the pointe we just located.
    while any(isnan(idx))
        counter = counter+1; 
        % find closest coordinate to (x(i),y(i)) that isn't already taken
        [~, idx(counter+1)] = min(pdist2(xyTemp,[col(idx(counter)), row(idx(counter))])); 
        % Eliminate that from the pool
        xyTemp(idx(counter+1),:) = NaN; 
    end

%---raw centerline in pixels
pxCenterline = [col(idx) row(idx)];


%---Find the extremities of the body on the outline---
% WARNING: the extremities of the centerline produced by BWskel do not land
% on the outline. As such, we need to extend the extremities of the
% centerline so that they eventually cross the centerline at the
% extremities. This process is iterative, but is very effective at finding
% the outline ends.
    % Extremity 1
        % find mean vector defining the direction of Extremity 1 in the x-y
        % plane. This is done by performing an average of the x and y vector
        % components of n segments spaning from the centerline tip to incrementally more distant points. This is done to find the general direction of the section while not being as sensitive to pixel shift. 
        nPt1 = 40; % number of points away from Extremity 1 in pixel. Change based on the resolution of the original video. The higher the resolution, the higher the number. WARNING: this assumes that the shape of the centerline along these n points is relatively straight.
        rawVecExt1 = zeros(nPt1,2); % store the x and y vectors for the direction of Extremity 1
        normRawVecExt1 = rawVecExt1;
            for ii = 1:nPt1
                rawVecExt1(ii,1) = pxCenterline(1,1)-pxCenterline(ii+1,1); % store the raw x vector component of the increasingly longer segment between Extremity 1 and nPt1
                rawVecExt1(ii,2) = pxCenterline(1,2)-pxCenterline(ii+1,2); % store the raw y vector component of the increasingly longer segment between Extremity 1 and nPt1
            end
        normRawVecExt1(:,1) = rawVecExt1(:,1)./(abs(rawVecExt1(:,1))+abs(rawVecExt1(:,2))); % normalize x vector
        normRawVecExt1(:,2) = rawVecExt1(:,2)./(abs(rawVecExt1(:,1))+abs(rawVecExt1(:,2))); % normalize y vector

        meanVecExt1 = mean(normRawVecExt1,1); % mean direction inferred from nPt1 points, assuming they are relatively alined in a straight line

        % Iteratively determine when the prolongation of Extremity 1 crosses the outline; move 1 px at a time
        nExt1 = 0; % initiate to stop the while loop when the condition is met
        iterExt1 = 1; % Initiate iteration. This number will grow until reaching a point when the centerline intersects with the outline.
            while nExt1 == 0
            [x1,y1] = polyxpoly(outline(:,1),outline(:,2),[pxCenterline(1,1) pxCenterline(1,1)+iterExt1*meanVecExt1(1,1)],[pxCenterline(1,2) pxCenterline(1,2)+iterExt1*meanVecExt1(1,2)]); % coordinates of Extremity 1 on the outline
            iterExt1 = iterExt1+1; % continue iteration
            if isempty(x1)
            nExt1 = 0;
            else
             nExt1 = 1; % stops the loop
            end
            end

    % Extremity 2
        % find mean vector defining the direction of Extremity 2 in the x-y
        % plane. This is done by performing an average of the x and y vector
        % components of n segments spaning from the centerline tip to incrementally more distant points. This is done to find the general direction of the section while not being as sensitive to pixel shift. 
        nPt2 = 20; % number of points away from Extremity 1 in pixel. Change based on the resolution of the original video. The higher the resolution, the higher the number. WARNING: this assumes that the shape of the centerline along these n points is relatively straight.
        rawVecExt2 = zeros(nPt2,2); % store the x and y vectors for the direction of Extremity 1
        normRawVecExt2 = rawVecExt2;
            for ii = 1:nPt2
                rawVecExt2(ii,1) = pxCenterline(end,1)-pxCenterline(end-ii,1); % store the raw x vector component of the increasingly longer segment between Extremity 2 and nPt2
                rawVecExt2(ii,2) = pxCenterline(end,2)-pxCenterline(end-ii,2); % store the raw y vector component of the increasingly longer segment between Extremity 2 and nPt2
            end
        normRawVecExt2(:,1) = rawVecExt2(:,1)./(abs(rawVecExt2(:,1))+abs(rawVecExt2(:,2))); % normalize x vector
        normRawVecExt2(:,2) = rawVecExt2(:,2)./(abs(rawVecExt2(:,1))+abs(rawVecExt2(:,2))); % normalize y vector

        meanVecExt2 = mean(normRawVecExt2,1); % mean direction inferred from nPt2 points, assuming they are relatively alined in a straight line

        % Iteratively determine when the prolongation of Extremity 2 crosses the outline; move 1 px at a time
        nExt2 = 0; % initiate to stop the while loop when the condition is met
        iterExt2 = 1; % Initiate iteration. This number will grow until reaching a point when the centerline intersects with the outline.
            while nExt2 == 0
            [x2,y2] = polyxpoly(outline(:,1),outline(:,2),[pxCenterline(end,1) pxCenterline(end,1)+iterExt2*meanVecExt2(1,1)],[pxCenterline(end,2) pxCenterline(end,2)+iterExt2*meanVecExt2(1,2)]); % coordinates of Extremity 2 on the outline
            iterExt2 = iterExt2+1; % continue iteration
            if isempty(x2)
            nExt2 = 0;
            else
            nExt2 = 1; % stops the loop
            end
            end

%---Concatenate the original BWskel centerline and the coordinates of the extremities
pxCenterlineFull = [x1,y1;pxCenterline;x2,y2];

% %---Create centerline with evenly distributed points---
% until the final step in the master script
N = 1000; % number of interpolated points; Use a number that is easy to normalize to identify a particular point along the segment. If N = 100 and th center of mass is at 0.25BL, then the x-y indices of CM will be 25.
centerline = curvspace([pxCenterlineFull(:,1) pxCenterlineFull(:,2)],N);
% centerline = [pxCenterlineFull(:,1) pxCenterlineFull(:,2)];

if opts.checkCenterline == 1
hold on
imagesc([0,(size(BW,2)-1)],[0,(size(BW,1)-1)],BWtempDilate); colormap('gray') % plot BW mask
plot(outline(:,1), outline(:,2), '-r','LineWidth',1);

% Plot the centerline
plot(centerline(:,1), centerline(:,2), '-b','LineWidth',2);

% Identify the ends
plot(centerline(1,1),centerline(1,2),'.g','MarkerSize',15) % Extremity 1
plot(centerline(end,1),centerline(end,2),'.g','MarkerSize',15) % Extremity 2

% Plot the center of mass
% CM = 0.25; % change accordingly
% plot(centerline(CM*N,1),centerline(CM*N,2),'oc')

% Label the axes to remember the scale
xlabel('distance (px)')
ylabel('distance (px)')
axis([0 (size(BW,2)-1) 0 (size(BW,1)-1)])
axis equal
end

end