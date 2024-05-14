%% Compute centerline from outline and calculate heading
% Author: Nils Tack
% Compiled with Matlab 2022a
% Latest update: 04-23-2024
% Script used to compute the centerline of an object and calculate its
% heading. This is useful for the computation of thrust drag, which
% requires a heading to perform accurate computations.
% 
% Before running the script, be sure to have access to the BW masks and
% outlines.

%% Clear workspace
close all; clearvars

%% Set paths
path.main = 'G:\Matlab scripts\computePlotPressureforces\Data'; % set path to 'Data' folder

% Other paths
path.masksBW = fullfile(path.main,'masksBW');       % BW masks produced from original images
path.outlines = fullfile(path.main,'outlines');     % outlines with x and y coordinates produced from BW masks (.csv)

% Set directories
D.masksBW = dir(fullfile(path.masksBW,'*.tif'));            % set BW masks files directory
D.outlines = dir(fullfile(path.outlines,'*.csv'));          % set outline files directory

%% -----LOAD THE SCALE-----
% Use the scale used to make the outlines and perform PIV processing in Davis
opts.scale = 15.23707*1000; % was measured in px/mm with ImageJ but needs to be to px/m for the rest of the calculations

%% -----LOAD MASKS-----
% For convenience, either the BW masks or the outlines can be used to
% compute the centerlines. Some people clean up the BW mask to opbtain
% smooth masks while other work on the outline instead.
% No conversion of the BW masks is needed, but if yousing the outlines,
% then the outlines are converted into BW masks so the skeletonize function
% performing the centerline computation can work
opts.BWorOutline = 2; % Two options: 1 = load the BWmasks images; 2 = load the outlines in xy coordinates and convert to BWmasks for centerline computation

if opts.BWorOutline == 1 % load BW mask
Data.BW = importMasks(D.masksBW); % import the BW masks
% flip the masks up-down to flip the orizontal axis so that the origin is
% at the bottom left of the frame

elseif opts.BWorOutline == 2 % convert outlines to BW masks
Data.BW = outline2masks(D,opts);
end

%% -----COMPUTE THE CENTERLINES-----
% Options
opts.checkCenterline = 1; % option to check the centerlines; 0 = disabled, 1 = enabled

%% Test one centerline
i = 1;
if opts.checkCenterline == 1
figure; hold on
end
testCenterline = makeCenterline(Data.BW(:,:,i),opts);
title(sprintf('Centerline %i',i))

%% Compute all the preliminary centerline
% the anterior extremity needs to be fixed because the centerline
% computation occasionally permutes the extremities. This is delt with in
% the next step.
if opts.checkCenterline == 1
fig = figure; hold on
pause(0.5) % give time for figure to load
end

for i = 1:length(D.masksBW)
Data.prelimCenterlines(:,:,i) = makeCenterline(Data.BW(:,:,i),opts);
    if opts.checkCenterline == 1
        title(sprintf('Centerline %i',i))
        pause(0.1) % give time for plot to update
    clf % clear figure so we don't overlap sevral together
        if i == length(D.masksBW)
            close(fig)
        end
    end
end

%% Finalize the certerlines and extract the heading of the fish
% 1) Place the coordinates of the snout of the fish as the first point in the
% centerline matrix
% 2) Compute the swimming direction of the fish assuming the head
% determines the heading. The heading is relative to true north in the
% image (+y direction).
Data = cleanCenterlines(Data,opts);

%% Convert the centerlines from px to meter
Data.centerlinesMeter = Data.centerlines/opts.scale;

%% Export centerlines and heading
exportCenterlines(path,Data)

%% Example plot showing everything
i = 1; % test frame
figure; hold on
imagesc([0,(size(Data.BW(:,:,i),2)-1)/opts.scale],[0,(size(Data.BW(:,:,i),1)-1)/opts.scale],Data.BW(:,:,i)); colormap('gray') % plot mask using BW image directly
plot(Data.centerlinesMeter(:,1,i),Data.centerlinesMeter(:,2,i),'-r','LineWidth',2) % plot the centerline
plot(Data.centerlinesMeter(1,1,i),Data.centerlinesMeter(1,2,i),'.r','MarkerSize',20) % plot the snout

% Plot the center of mass
% The number of points along the centernile is set to 1000 by default (see
% 'makeCenterline', line 123).

CM = 0.30; % distance of the CM from the snout as a fraction of body length, in %
CMindex = [Data.centerlinesMeter(CM*1000,1,i),Data.centerlinesMeter(CM*1000,2,i)];
plot(CMindex(:,1),CMindex(:,2),'.k','MarkerSize',20)

title(sprintf('Centerline %i',i))
xlabel('distance (m)')
ylabel('distance (m)')
% axis([0 (size(Data.BW(:,:,i),2)-1) 0 (size(Data.BW(:,:,i),1)-1)])
axis equal