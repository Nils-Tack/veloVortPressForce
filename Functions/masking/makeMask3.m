function BW5 = makeMask3(I, opts)
if opts.inverseImage % option to inverse the image so that the object to be masked appears white
   I = imcomplement(I); % inverse the image;
end

%------CONTRAST ENHANCEMENT-----
% Adjust contrast and gamma correction on original image
I2 = imadjust(I,[0.40 1],[],0.9); % change values accordingly. (Image, [low_in high_in], [], gamma); If gamma is less than 1, then imadjust weights the mapping toward higher (brighter) output values.

%-----PERFORM BINARIZATION-----
BWopt = 1;          % Change option according to preferences. Binarization option: 0 = normal binarization, 1 = adaptive binarization.
BWthreshold = 0.45;  % Change threshold value for binarization based on the quality of the review. Number in the range [0, 1]. A high sensitivity value leads to thresholding more pixels as foreground, at the risk of including some background pixels.

if BWopt == 0
   BW2 = imbinarize(I2, BWthreshold);
else
   BW2 = imbinarize(I2,'adaptive','Sensitivity',BWthreshold); % Sensitivity factor for adaptive thresholding; higher = more details, lower = less details
end

%-----MANUALLY MASK HORIZONTAL AND VERTICAL SECTIONS OF THE IMAGE-----
% This may need to be performed to remove unwanted bright objects in the
% field of view that would throw off the detection of the object of
% interest in the next step.
Maskopt = 0; % Change according to whether you want to eliminate areas in the rame from analysis. Options: 0 = do not add additional masking; 1 = add additional vertical masking; 2 = add horizontal masking; 3 = add both horizontal and vertical masking

if Maskopt == 1 % vertical masking
   BW2(:,1:600) = 0; % left; change indices as needed
   BW2(:,size(BW2,1)-400:size(BW2,1)) = 0; % right; change indices as needed

elseif Maskopt == 2  % horizontal masking
   BW2(1:300,:) = 0; % top; change indices as needed
   BW2(size(BW2,1)-200:size(BW2,1),:) = 0; % bottom; change indices as needed

elseif Maskopt == 3  % both horizontal and vertical masking
    BW2(:,1:600) = 0; % left; change indices as needed
    BW2(:,size(BW2,1)-400:size(BW2,1)) = 0; % right; change indices as needed
    BW2(1:300,:) = 0; % top; change indices as needed
    BW2(size(BW2,1)-200:size(BW2,1),:) = 0; % bottom; change indices as needed
end

%-----CLEAN THE MASK-----
% Perform morphological operations to eliminate unwanted features

% 1) Extract the largest blob(s) only
BlobOpt = 2; % Option to save the largest object only, or all the objects above a certain size; Options: 1 = largest object only, 2 = Any object larger than a threshold.

% Check that the frame contains at least one non zero element - use the sum
% of all pixel.
if sum(BW2(:))>0 % necessary to declare in case there are no objects in the field of view. This will enabe the production of a completely white mask that won't disrupt subsequent plots.
props = regionprops(BW2, 'Area'); % Find areas of each blob and put into a structure.
allAreas = [props.Area];
    if BlobOpt == 1
       BW3 = bwareafilt(BW2, 1); %  keeps the n largest objects; change based on the number of masks needed in one image
    elseif BlobOpt == 2
       blobSize = 500; % minimum threshold; change as necessary
       if max(allAreas)>=blobSize
       BW3 = bwareafilt(BW2, [blobSize,max(allAreas)]); % select all the blob ranging from the largest to a specific, smaller threshold
       end
    end
if max(allAreas)>=blobSize
% 2) Fill holes
BW3 = ~bwareaopen(~BW3, 10); % awesome to only fill small holes. Perfect for complicated shapes and systems like multiple legs/appendages; change value as needed.
BW3 = imfill(BW3, 'holes'); % Closes larger holes

% 3) Smooth edges (through aggressive, temporary dillation)
% Dilate the Image to preserve small details (like the tip of the tail)
se90 = strel('line',10,90);
se0 = strel('line',10,0);

% 4) Thicken the mask
BW4 = bwmorph(BW3,'thicken',5); % was 10

% 5) Blur the image
windowWidth = 15; % was 25
BWblurr = conv2(double(BW4), ones(windowWidth)/windowWidth^2, 'same');

% 6) Threshold again.
BWblurr = BWblurr > 0.5;

% 7) Shrink smooth mask to original size
BW5 = bwmorph(BWblurr,'shrink',2.5); % 1/2 of original dilation since blur was cut to 1/2 already
else
BW5 = BW2*0;
end


% Plot combination of original, blurred, and final mask images
if opts.plotComparison
f = figure;
f.Position(1:2) = [0 0];
subplot(2,2,1)
imshow(I)
title('Original image')
subplot(2,2,2)
imshow(I2)
title('Gamma correction')
subplot(2,2,3)
imshow(BW2)
title('Adaptive binarization + removing edges')
subplot (2,2,4)
imshow(BW5);
title('Smoothed mask')
end

end