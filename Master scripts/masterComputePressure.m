%% Make masks, outlines, and compute pressure
% Author: Nils Tack
% Compiled with Matlab 2023b
% Latest update: 04-29-2024
% Script used to prepare PIV data files, create masks/outlines from images,compute pressure fields, thrust
% and drag forces acting on the body.

% Please refer to the following publication for details about pressure field 
% calculations from 2D velocity fields:
% Dabiri, J. O., Bose, S., Gemmell, B. J., Colin, S. P. and 
% Costello, J. H. (2014). An algorithm to estimate unsteady and quasi-steady 
% pressure fields from velocity field measurements. J. Exp. Biol. 217, 331–336.

% Before running the script, copy raw PIV data text files into the 'PIVoriginal' folder
% found in 'Data'
% To make masks, also load corresponding image sequence into the 'images' folder
% found in 'Data'. Do not worry about frame increment yet. Copy all images based on
% start and end frames corresponding to the file number of first and last
% PIV text files.

% Required file format: image files (.jpg or .tif), PIV data (.txt)
% WARNING: .png images can cause issues with the colormap.
%% Clear workspace
close all; clearvars

%% Set paths
path.main = 'G:\Matlab scripts\computePlotPressureforces\Data'; % set path to 'Data' folder

% Other paths
path.PIVraw = fullfile(path.main,'PIVoriginal');    % raw PIV data
path.images = fullfile(path.main,'images');         % image sequence
path.PIVclean = fullfile(path.main,'PIVclean');     % clean velocity data without headers (convert data for use by Queen2)
path.masksBW = fullfile(path.main,'masksBW');       % BW masks produced from original images
path.outlines = fullfile(path.main,'outlines');     % outlines with x and y coordinates produced from BW masks (.csv)
path.pressure = fullfile(path.main,'pressure');     % pressure data (.csv)
path.forces = fullfile(path.main,'forces');         % force data

% Set file type for image sequence (.jpg or .tif or .png)
opts.imFileType = '.jpg'; % !!! warning: png cause issues with colormap overlay.

% Set directories (only for existing files)
D.images = dir(fullfile(path.images,['*',opts.imFileType])); % make image directory
D.PIVraw = dir(fullfile(path.PIVraw,'*.txt')); % set raw PIV files directory

%% -----CONVERT RAW VELOCITY FILES TO .CSV-----
txt2csv(D,path)

%% -----MASKING-----
%% Options for masks
% select the preferences for a given project. 1 enables the function, 0
% disables it.
opts.inverseImage    = 0; % inverse the image for masking (negative). This is useful when working with laser pIV and bright field PIV. Masking assumes that the obesct to mask is white (logical  = 1) and the background black (logical = 0)
opts.checkMasks     = 1; % check the quality of the masks as they get exported. Useful to determine whether the masks will need some manual processing.
opts.maskImage      = 0; % Point-click masking tool to pre-process the image and black out unwanted areas of the frame that may throw off the masking script. Useful to remove the sides/walls of a tank. This function uses a derived version of ginput. WARNING! this function is useful but not recommended since it will not automatically save the coordinates of the points making up the masked portion of the image.
opts.plotComparison = 1; % plots mask-making steps for troubleshooting during the initial mask testing phase (muted automatically for mask export)

%% Mask fish - test frame
% Read a test image
f = 10; % test frame number
I = importFrame(D.images,f); % import one frame
filename = fileID(D.images,f); % extract file name for use as figure title

if opts.maskImage % Click and drag option to select only the area of interest. This automatically masks unnecessary parts of the image outside the selection
    [I,x_ginput,y_ginput] = maskMasked(I); % select two oposing corners of a rectangle including the area of interest. This function returns the image with black masked out areas that can be used directly with makeMask3
end

% Open 'makeMask3' function to change parameters if necessary
% Parameters include masking unwanted areas like tank/flume edges or
% reflections, contrast and gamma corrections, smoothing options.
I4 = makeMask3(I, opts);

%% Make and export all the masks
opts.plotComparison = 0; % plots mask-making steps (muted automatically for mask export)

% Load all the images into a 3D matrix
Data.images = importFrame(D.images); % import all the frames; this step can take a few seconds depending on the size of your image stack
Data.filename = fileID(D.images); % extract all the file name for use as figure title

if opts.checkMasks % check individual masks as they get exported
    figTemp = figure;
end

fprintf('Exporting masks...');
for i = 1:length(D.images)
    I = Data.images(:,:,i); % load corresponding image;

    % Masking unnecessary parts of the image to isolate the subject based
    % on selection above test section.
    if opts.maskImage
       I([1:y_ginput(1),y_ginput(2):end], :) = 0;
       I(:, [1:x_ginput(1),x_ginput(2):end]) = 0;
    end

% Make a binary mask
I4 = makeMask3(I,opts);
    
% Option to check masks during export
    if opts.checkMasks
       imshow(I4)
       title(Data.filename{i}) % WARNING: this will display the actual file name of the corresponding image. However, the masks will be stored with names atrting at 00001!
       pause(0.05);
       clf
    end
        
% Export mask
filenameBW = sprintf('BW_%05g',i);
imwrite(I4,fullfile(path.masksBW,[filenameBW,'.tif'])) % exports as tif     
end

if opts.checkMasks % check individual masks as they get exported
close(figTemp)
end

fprintf('done\n');

%% -----GENERATE OUTLINE-----
% Set scale to make the outlines (same as that used for PIV)
opts.scale = 15.8372*1000; % was measured in px/mm with ImageJ but needs to be to px/m for the rest of the calculations

% Options
opts.multiOutlines = 0; % 0 = uses only the biggest blob created during the BW mask export step; 1 = option to enable the storage of multiple outlines in the same file

% set BW mask directory
D.masksBW = dir(fullfile(path.masksBW,'*.tif')); 

% Import all masks
Data.masksBW = importMasks(D.masksBW); % import all the masks


%% Test one outline
figure; hold on % initiate figure
outline = makeOutlines(Data.masksBW(:,:,f),opts); % make a test mask
title(sprintf('Outline %05g',f));

%% Export all the outlines
figure; hold on % initiate figure

fprintf('Exporting outlines...');
for i = 1:size(Data.masksBW,3)
outline = makeOutlines(Data.masksBW(:,:,i),opts); % make a test mask
title(sprintf('Outline %05g',i));

% export individual outline to a .csv file
filenameOutline = sprintf('iface_%05g', i); % Number sequentially
writematrix(outline,fullfile(path.outlines,[filenameOutline,'.csv']))
progressCount(i,size(Data.masksBW,3)); % display export progress

pause(0.01)
hold off
clf
end
close
fprintf('done\n');

%% -----OPTIONAL - REPAIR OUTLINES MANUALLY IF NECESSARY-----
% Set raw outline directory
D.outlines = dir(fullfile(path.outlines,'*.csv')); 

% Repair outlines manually
repairOutlines(D,[1 length(D.outlines)],opts) % change the first and last frame accordingly [first last].

%% Export BW masks with all the clean outlines (OPTIONAL)
exportCleanMasks(path,D,size(Data.masksBW),opts)

%% -----CALCULATE PRESSURE FIELDS-----
% Open parameters used by queen2 to run computation
% Change the required fields in the struct element
open('parameters2')

%% Run queen2
% Updated colormap for live pressure fields plot to enhance readability
% (see custom colormap lines 421 & 1090 for details)
queen2

%% Move pressure files from the PIVclean folder to pressure folder
D.preliminaryPress = dir(fullfile(path.PIVclean,'press*.csv'));
movePress(D,path)

