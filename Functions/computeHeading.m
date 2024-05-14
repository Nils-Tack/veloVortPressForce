function heading_degrees = computeHeading(pxCenterline)
% Calculate the heading (swimming direction) of the fish in degrees relative to the +y axis (similar to true north when calculating the heading)

nPt1 = 150; % number of points away from the tip of the snout. Change based on the resolution of the original video. The higher the resolution, the higher the number. WARNING: this assumes that the shape of the centerline along these n points is relatively straight. This should be the case for the head of a fish that is completely stiff.
rawVecExt1 = zeros(nPt1,2); % store the x and y vectors for the direction of Extremity 1
normRawVecExt1 = rawVecExt1;
    for ii = 1:nPt1
        rawVecExt1(ii,1) = pxCenterline(1,1)-pxCenterline(ii+1,1); % store the raw x vector component of the increasingly longer segment between Extremity 1 and nPt1
        rawVecExt1(ii,2) = pxCenterline(1,2)-pxCenterline(ii+1,2); % store the raw y vector component of the increasingly longer segment between Extremity 1 and nPt1
    end
normRawVecExt1(:,1) = rawVecExt1(:,1)./(abs(rawVecExt1(:,1))+abs(rawVecExt1(:,2))); % normalize x vector
normRawVecExt1(:,2) = rawVecExt1(:,2)./(abs(rawVecExt1(:,1))+abs(rawVecExt1(:,2))); % normalize y vector

meanVecExt1 = mean(normRawVecExt1,1); % mean direction inferred from nPt1 points, assuming they are relatively alined in a straight line

% Define the direction vector components (x, y)
direction_vector = [meanVecExt1(2), meanVecExt1(1)]; % WARNING: load as [y,x]

% Calculate the heading angle in radians
heading_radians = atan2(direction_vector(2), direction_vector(1));

% Convert radians to degrees
heading_degrees = rad2deg(heading_radians);

% Ensure the heading angle is in the range [0, 360)
if heading_degrees < 0
    heading_degrees = heading_degrees + 360;
end
% 
% disp(['Heading angle in degrees: ', num2str(heading_degrees)]);


end
