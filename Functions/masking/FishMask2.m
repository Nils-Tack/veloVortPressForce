function [Mask] = FishMask2(Frames)

 [a b c] = size(Frames);
parfor e = 1:c

    Io = Frames(:,:,e);     % Original Image
    s = min(Io,[],'all');       % Lowest intensity in the image
    Si = Io == s;           % Binary image with the lowest intensity pixels (i.e., animal detection)

    for k = 1:50                % This for loop adds pixels with intensities up to +50 from the darkest point;
        s = min(min(Io))+k;     % It is still dark and makes the seed a bit larger.
        So = Io == s;
        Si = Si+So;
    end

    S = bwareaopen(Si,1000);        % Final seed image for regiongrow function. 
    % bwareaopen get rid of areas that are two small (probably not inside
    % tha animal). Modify this threshold as needed. 


% -------------------------------------------------------------------------
% Here ROI is a region of interest. By combining S with I1 and agressively
% dilating I'm getting a region that for sure has the compleate animal
% inside. 

    I1 = edge(Io,'Canny',[0.1 0.7]);

    ROI = imbinarize(S+I1);
    se = strel('disk',5,8);
    ROI = imclose(ROI,se);
    
    aa = bwarea(ROI);
    ROI = bwareaopen(ROI,round(0.01*aa));      

    se2 = strel('disk',50,8);
    ROI = imdilate(ROI,se2); 

% Then I apply a gaussian filter that gets rid of the fluid particles in
% this ROI, followed by a contrast enhacement. All these steps help for a
% more accurate animal detection later on. 
    f = @(Io) imgaussfilt(Io,20,'FilterSize',3);
    I1 = roifilt2(Io,ROI,f);

    L = stretchlim(Io);
    f1 = @(I1) imadjust(I1,L,[],0.5);
    I1 = roifilt2(I1,ROI,f1);


% Finally I use the regiongrow to find the animal. 

    To = 50; 
    f2 = @(I1) regiongrow(I1,S,To);
    I2 = roifilt2(I1,ROI,f2);
    I2 = I2.*ROI;
    I2 = imdilate(I2,se);
    A = bwarea(I2);
    I2 = ~bwareaopen(~I2,round(0.001*A));

 % 
 figure(1);
    subplot(2,2,1)
    imshow(Io)
    subplot(2,2,2)
    imshow(S)
    subplot(2,2,3)
    imshow(I1)
    subplot(2,2,4)
    imshow(I2)

    Mask(:,:,e) = I2;
end

end