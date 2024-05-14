%% Plot velocity, vorticity, and pressure fields
% Author: Nils Tack
% Compiled with Matlab 2022a
% Latest update: 04-18-2024
% Script used to prepare clean plots of the velocity, vorticity, and
% pressure fields.
% 
% Before running the script, be sure to have access to the PIVclean
% velocity files (.csv format), and the pressure files (.csv). Also have
% the BW or outlines ready to overlay the masks
% 
% Required file format: image files (.jpg or .tif ), PIV and pressure data (.csv)
% WARNING: .png images can cause issues with the colormap.

%% Clear workspace
close all; clearvars

%% Set paths
path.main = 'G:\Matlab scripts\computePlotPressureforces\Data'; % set path to 'Data' folder

% Other paths
path.images = fullfile(path.main,'images');         % image sequence
path.PIVclean = fullfile(path.main,'PIVclean');     % clean velocity data without headers (convert data for use by Queen2)
path.masksBW = fullfile(path.main,'masksBW');       % BW masks produced from original images
path.outlines = fullfile(path.main,'outlines');     % outlines with x and y coordinates produced from BW masks (.csv)
path.pressure = fullfile(path.main,'pressure');     % pressure data (.csv)


% Set file type for image sequence (.jpg or .tif or .png)
opts.imFileType = '.jpg'; % !!! warning: png cause issues with colormap overlay.

% Set directories
D.images = dir(fullfile(path.images,['*',opts.imFileType]));     % make image directory
D.PIVclean = dir(fullfile(path.PIVclean,'*.csv'));          % set raw PIV files directory
D.masksBW = dir(fullfile(path.masksBW,'*.tif'));            % set BW masks files directory
D.outlines = dir(fullfile(path.outlines,'*.csv'));          % set outline files directory
D.pressure = dir(fullfile(path.pressure,'*.csv'));  % set pressure files directory

%% -----LOAD THE SCALE-----
% Use the scale used to make the outlines (same as that used for PIV processing in Davis)
opts.scale = 15.8372*1000; % was measured in px/mm with ImageJ but needs to be to px/m for the rest of the calculations

%% -----VELOCITY FIELDS-----
%% Import and store all the velocity fields
veloData = importVeloFields(D); % structure containing meshgrid X and Y for velocity, vorticity, and pressure plots; The structure also stores the u and v data as well as the absolute speed for each X,Y point.

%% Convert sections falling outside the computation domain to NaN
% Useful if a geometric mask was applied in DaVis. It is necessary to tuen 0s to NaN to plot clean vorticity fields.
veloData = veloZeros2NaN(veloData);
% KEEP IN MIND: Davis will never return a strict 0 velocity value, thus we
% can safely assume that any of the true 0 u and v values fell outside the
% domain of interest in Davis, and can thus be considered NaN.

%% Plot mean velocity field and subtract bulk flow if necessary
% Necessary to measure the bulk flow when using a flume.
% Subtracting background flow from u or v shows vortices. If no bulk flow, no subtraction is performed

veloData.meanU = mean(veloData.uVelo,3); % mean u performed across the 3D matrix
veloData.meanV = mean(veloData.vVelo,3); % mean v performed across the 3D matrix

% Plot mean velocity field
% Plot the time-averaged u or v component to measure the mean bulk flow
opts.flowDirection = 2;         % 1 = u (x direction); 2 = v (y direction); 3 = no flow
opts.flowVeloALreadyKnown = 0;  % 0 = the velocity of the flow is not already known; 1 = the velocity of the flow in u or v is already known and replaces the need to measure it.
opts.fixedVelo = 0.16;  % Fixed volocity to subtract if already known. Applicable to opts.flowDirection = 1 or 2.

veloData = plotMeanVelo(veloData,opts); % keep in mind that subtraction of the bulk flow (if applicable) is performed here.

%% -----PLOT VELOCITY AND VORTICITY FIELDS-----
% Options
opts.plotVectors = 1;       % Option to overlay the velocity vectors 0 = no vector; 1 = vector fields (preferred)
opts.plotMask = 0;          % Plot masks (instead of leaving the vorticity field empty) 0 = no overlaid mask (preferred to see the animal); 1 = plain mask (black or white)
opts.tightPlot = 0;         % Option to plot specific section of the velocity/vorticity fields (good for making clean figures, but the cropping causes issues with the exported .svg files in Illustrator. I recommend keeping 0
opts.outline = 1;           % Option to use the outlines rather than the BW masks to mask the velocity data. 0 = do not utilize either BW masks or outlines; 1 =  use existing outlines (preferred), 2 = do not utilize outlines and use the raw BW masks instead (could be useful sometimes
opts.exportPlot = 1;        % options to chose the figure export format. Three options: 1) Export only jpg; 2) export only .svg; 3) export both .jpg and .svg.

% Plotting options
opts.vortMag = 20;          % Option to set the colorbar min and max vorticity magnitude
opts.Resolution = 2;        % Option to change the resolution in two different ways. 1 = set the increment between vectors by defining the actual physical spacing in meter; 2) Increase the resolution by a factor of n. The higher the resolution, the longer it will tekae to interpolate the fields.
opts.smoothness = 5;        % Spacing between vectors; use value corresponding to the above 'Resolution' option
opts.vecDensity = 1;        % Parameter to plot 1/n vector. 1 plots all the vectors. The higher the number, the sparser the velocity field. Use integer only.

warning('off','MATLAB:polyshape:repairedBySimplify');   % disable unecessary polyshape warning
warning('off','MATLAB:polyshape:boundary3Points');      % disable unecessary polyshape warning

%% Run a test frame
f = 20; % test frame ID
plotVeloVort(veloData,[f f],opts,path,D) % velocity data, [frames for export, first and last], options, paths, directories

%% Export all the velocity and vorticity fields
plotVeloVort(veloData,[1 size(veloData.uVelo,3)],opts,path,D)

%% -----PLOT PESSURE-----
% Options
opts.plotMask = 0;          % Plot masks (instead of leaving the vorticity field empty) 0 = no overlaid mask (preferred to see the animal); 1 = plain mask (black or white)
opts.tightPlot = 0;         % Option to plot specific section of the velocity/vorticity fields (good for making clean figures, but the cropping causes issues with the exported .svg files in Illustrator. I recommend keeping 0
opts.outline = 1;           % Option to use the outlines rather than the BW masks to mask the velocity data. 0 = do not utilize either BW masks or outlines; 1 =  use existing outlines (preferred), 2 = do not utilize outlines and use the raw BW masks instead (could be useful sometimes
opts.exportPlot = 1;        % options to chose the figure export format. Three options: 1) Export only jpg; 2) export only .svg; 3) export both .jpg and .svg. (extremely slow)

% Plotting options
opts.pressMag = 15;         % Option to set the colorbar min and max pressure magnitude. The value you input is the min and max of the colorbar.
opts.Resolution = 2;        % Option to change the resolution of the vorticity fields in two different ways. 1 = set the increment between vectors by defining the actual physical spacing in meter; 2) Increase the resolution by a factor of n. The higher the resolution, the longer it will tekae to interpolate the fields.
opts.smoothness = 5;        % Spacing between vectors; use value corresponding to the above 'Resolution' option

warning('off','MATLAB:polyshape:repairedBySimplify');   % disable unecessary polyshape warning
warning('off','MATLAB:polyshape:boundary3Points');      % disable unecessary polyshape warning

%% Import and store all the velocity fields
pressData = importPressFields(D); % structure containing meshgrid X and Y for velocity, vorticity, and pressure plots; The structure also stores the u and v data as well as the absolute speed for each X,Y point.

%% Run a test frame
f = 20; % test frame ID
plotPress(pressData,[f f],opts,path,D) % velocity data, [frames for export, first and last], options, paths, directories

%% Export all the pressure fields
plotPress(pressData,[1 size(pressData.P,3)],opts,path,D)