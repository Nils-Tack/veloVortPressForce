%% Compute and plot forces from pressure fields
% Author: Nils Tack
% Compiled with Matlab 2023b
% Latest update: 04-29-2024
% Script used to compute and plot force vectors thrust, drag, and lateral
% forces) using the heading of the object in the fluid.

% Please refer to the following publication for details about force
% calculations:
% Tack, N. B., Du Clos, K. T., & Gemmell, B. J. (2024). Fish can use coordinated fin motions to recapture their own vortex wake energy. Royal Society Open Science, 11(1), 231265.

% Before running the script, make sure to have the pressure fields and
% outlines of the object in the fluid. Also make sure to have the heading
% the the object, either as a pre-determined unchanging value, or as a single .csv
% file (n x 2 matrix) exported when acquiring the centerlines.

%% Clear workspace
close all; clearvars

%% Set paths
path.main = 'G:\Matlab scripts\computePlotPressureforces\Data'; % set path to 'Data' folder

% Other paths
path.images = fullfile(path.main,'images');             % image sequence
path.masksBW = fullfile(path.main,'masksBW');           % BW masks produced from original images
path.outlines = fullfile(path.main,'outlines');         % outlines with x and y coordinates produced from BW masks (.csv)
path.pressure = fullfile(path.main,'pressure');         % pressure data (.csv)
path.forces = fullfile(path.main,'forces');             % force data
path.centerlines = fullfile(path.main,'centerlines');   % centerlines data (also containing the swimming direction data)


% Set file type for image sequence (.jpg or .tif or .png)
opts.imFileType = '.jpg'; % !!! warning: png cause issues with colormap overlay.

% Set directories (only for existing files)
D.images = dir(fullfile(path.images,['*',opts.imFileType]));    % set image directory
D.masksBW = dir(fullfile(path.masksBW,'*.tif'));                % set BW masks files directory
D.outlines = dir(fullfile(path.outlines,'*.csv'));              % set outline files directory
D.pressure = dir(fullfile(path.pressure,'*.csv'));              % set pressure files directory
D.centerlines = dir(fullfile(path.centerlines,'*.csv'));        % set centerlines directory
D.heading = dir(fullfile(path.centerlines,'heading','*.csv'));  % set heading matrix directory

%% Options
opts.scale = 15.8372*1000;          % was measured in px/mm with ImageJ but needs to be to px/m for the rest of the calculations
opts.pressureDeltaT  = 0.0005;      % Time interval between pressure fields in second
opts.heading = 2;                   % Object heading (relative to +y, in degrees). 1 = fixed heading; 2 = changing heading (required matrix of length D.outlines).
opts.fixedHeading = 180;            % Fixed heading of the object relative to +y. Used only if opts.heading = 1; otherwise ignored by the force calculation function
opts.swimmingSpeed   = 1;           % Average swimming speed in meter per second of the animal assuming no acceleration (this can be updated inthe future by measuring the displacement of the center of mass from one frame to the next when exporting the heading)


% Plotting options
opts.npoints         = 300;         % Number of points along the body to trace force vectors (arbitrary)
opts.forceScale      = 2;           % Scale factor for plotting force vectors (for plotting purposes only, does not affect computation of forces)
opts.increment       = 1;           % Force calculation increment (from one frame to the next); generally stays at 1
opts.useRawImage     = 1;           % Option to plot raw image or mask (from outline); 0 = uses mask only (from outline), 1 = uses the raw image and overlays force vectors
opts.exportFigures   = 1;           % Option to export figures in the 'forces' folder
opts.exportPlotFormat = 1;          % options to chose the figure export format. Three options: 1) Export only jpg; 2) export only .svg; 3) export both .jpg and .svg.

%% Test one frame
f = 1;
% options for plot: 'axial', 'lateral', 'both', 'all'
[summaryDataForceTest] = forcesFromPressure(D,path,opts,'all',[f f]);

% Plot the force legend for reference
plotForceLegend
%% Export forces for all the frames
% options for plot: 'axial', 'lateral', 'both', 'all'
[summaryDataForceTest] = forcesFromPressure(D,path,opts,'axial',[1 length(D.pressure)]);
