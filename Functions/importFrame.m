function I = importFrame(D_images,f)
channelID = 1; % RGB channel of the image; check camera specs for which layer contains the most information. Most sensors are more sensitive to green, but for PIV, this will depend on the wavelength of the laser

if nargin > 1
% Import image data for one image # f
Itemp = importdata(quickfilepath(D_images(f)));
I = Itemp(:,:,channelID); % select only one of the three layers of the image;

else % import the full stack of image and store in a 3D array
    for i = 1:length(D_images)
        Itemp = importdata(quickfilepath(D_images(i)));
        channelID = 1;
        I(:,:,i) = Itemp(:,:,channelID);
    end

end
end