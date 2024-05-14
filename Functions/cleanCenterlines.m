function Data = cleanCenterlines(Data,opts)
% Parameters to control the size of the circular search area around the snout
r = 100; % radius of the seach area in px
n = 25; % number of points along the circumference of the search area(keep low since we just do a basic search here)

% Make temporary placeholder for the coordinates of the center of the
% search area
searchCenter = [0 0];

figure; 

for i = 1:size(Data.BW,3)
    hold on
    % Label the axes to remember the scale
    title(sprintf('Centerline %i',i))
    xlabel('distance (px)')
    ylabel('distance (px)')
    axis([0 (size(Data.BW(:,:,i),2)-1) 0 (size(Data.BW(:,:,i),1)-1)])
    axis equal

    imagesc([0,size(Data.BW(:,:,i),2)-1],[0,(size(Data.BW(:,:,i),1)-1)],Data.BW(:,:,i)); colormap('gray') % plot BW mask
    plot(Data.prelimCenterlines(:,1,i), Data.prelimCenterlines(:,2,i), '-b','LineWidth',2);  % Plot the first preliminary centerline

    if i < 3 % assign the search area to find the coordinates of the snout manually
        if i == 1
        % Identify the snout by selecting a region, use ginpuc to find a
        % rough area where to search. Only click once to define the center of the search area
        [searchCenter(1),searchCenter(2)] = ginputc(1,'Color','r'); % Find the coordinates of the center of the search area
        end

        % Make circular search area
        theta = (0:n-1)*(2*pi/n);
        xSearch = searchCenter(1) + r*cos(theta);
        ySearch = searchCenter(2) + r*sin(theta);
        plot(xSearch,ySearch,'.r','MarkerSize',10) % plot the search area

        % Find if the first point of the centerline is in the circle
        is_inside = inpolygon(Data.prelimCenterlines(1,1,i),Data.prelimCenterlines(1,2,i), xSearch, ySearch);
        
        if is_inside
            % The point is inside the search area
            Data.centerlines(:,:,i) = Data.prelimCenterlines(:,:,i); % keep the centerline the way it is  and store in D
        else
            % The point is outside the search area
            Data.centerlines(:,:,i) = flipud(Data.prelimCenterlines(:,:,i)); % flipup the centerline to make the snout the first point in the matrix and store in D
        end
        
        % plot the snout to ensure it is properly located
        plot(Data.centerlines(1,1,i),Data.centerlines(1,2,i),'.g','MarkerSize',15) % Snout

        % calculate the heading of the fish
        Data.heading(i,1) = computeHeading(Data.centerlines(:,:,i));
   
    else % assigns the search area based on the vector displacement of the snout in the previous frame. This requires good temporal and spatial resolution in the original video, but it is more reliable than simply using the coordinates of the previous frame to define the search area.
        % calculate the displacement vector of the snout from the two
        % previous frames to estimate the loaction of the search area in the
        % corresponding frame
        snoutDisp(1) = Data.centerlines(1,1,i-1)-Data.centerlines(1,1,i-2); % store the raw x vector component of displacement of the i-2 and i-1 frames to estimate the x coordinate of the search area in frame i
        snoutDisp(2) = Data.centerlines(1,2,i-1)-Data.centerlines(1,2,i-2); % store the raw y vector component of displacement of the i-2 and i-1 frames to estimate the y coordinate of the search area in frame i
        [searchCenter] = [Data.centerlines(1,1,i-1)+snoutDisp(1),Data.centerlines(1,2,i-1)+snoutDisp(2)];
        
        % Make circular search area
        theta = (0:n-1)*(2*pi/n);
        xSearch = searchCenter(1) + r*cos(theta);
        ySearch = searchCenter(2) + r*sin(theta);
        plot(xSearch,ySearch,'.r','MarkerSize',10) % plot the search area
        
        % Find if the first point of the centerline is in the circle
        is_inside = inpolygon(Data.prelimCenterlines(1,1,i),Data.prelimCenterlines(1,2,i), xSearch, ySearch);
        
        if is_inside
            % The point is inside the search area
            Data.centerlines(:,:,i) = Data.prelimCenterlines(:,:,i); % keep the centerline the way it is  and store in D

        else % we need to confirm whether the snout of the inverted profile fits within the search area or whether the search area was just misplaced due to to fast of a displacmenet between frames
            tempCenterline = flipud(Data.prelimCenterlines(:,:,i));
            is_inside_temp = inpolygon(tempCenterline(1,1),tempCenterline(1,2), xSearch, ySearch);

            if is_inside_temp
               Data.centerlines(:,:,i) = tempCenterline; % flipup the centerline to make the snout the first point in the matrix and store in D
            else
               [searchCenter(1),searchCenter(2)] = ginputc(1,'Color','r'); % Re-identify the coordinates of the center of the search area
                % Make circular search area
                theta = (0:n-1)*(2*pi/n);
                xSearch = searchCenter(1) + r*cos(theta);
                ySearch = searchCenter(2) + r*sin(theta);
                plot(xSearch,ySearch,'.r','MarkerSize',10) % plot the search area
                is_inside_temp2 = inpolygon(Data.prelimCenterlines(1,1,i),Data.prelimCenterlines(1,2,i), xSearch, ySearch);
                
               if is_inside_temp2
                    Data.centerlines(:,:,i) = Data.prelimCenterlines(:,:,i); % keep the centerline the way it is  and store in D
               else
                    Data.centerlines(:,:,i) = flipud(Data.prelimCenterlines(:,:,i)); % flipup the centerline to make the snout the first point in the matrix and store in D
               end

            end
            
        end

        % plot the snout to ensure it is properly located
        plot(Data.centerlines(1,1,i),Data.centerlines(1,2,i),'.g','MarkerSize',15) % Snout

        % calculate the heading of the fish
        Data.heading(i,1) = computeHeading(Data.centerlines(:,:,i));
    end

% % Label the axes to remember the scale
% xlabel('distance (px)')
% ylabel('distance (px)')
% axis([0 (size(Data.BW(:,:,i),2)-1) 0 (size(Data.BW(:,:,i),1)-1)])
% axis equal
pause(0.1)
hold off
clf

end
end
