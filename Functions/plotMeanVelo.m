function veloData = plotMeanVelo(veloData,opts)

if opts.flowVeloALreadyKnown == 0
figure; hold on
axis equal
vecScale = 0.1; % vector length scale
maxScale = 0.2; % max velocity magnitude
end

% Trace velocity profile across the field of view (run only if the flow speed is not known)
npoints = 200; % number of points along velocity profile to be extracted (for interpolation)

if opts.flowDirection == 1
   if opts.flowVeloALreadyKnown == 0
   quiverColor(veloData.X,veloData.Y,veloData.meanU,veloData.meanV*0,maxScale,vecScale);
   veloData.profileVelo = quickprofile(veloData.X,veloData.Y,veloData.meanU,npoints); % profile_u is the u velocity in the x direction (measured along y)
   % Subtract bulk flow from u
   veloData.uVelo = veloData.uVelo-mean(veloData.profileVelo);
   fprintf('Mean u flow speed: %f m.s-1\n',mean(veloData.profileVelo))
   elseif opts.flowVeloALreadyKnown == 1
   veloData.uVelo = veloData.uVelo-opts.fixedVelo;
   end
       

elseif opts.flowDirection == 2
    if opts.flowVeloALreadyKnown == 0
    quiverColor(veloData.X,veloData.Y,veloData.meanU*0,veloData.meanV,maxScale,vecScale);
    veloData.profileVelo = quickprofile(veloData.X,veloData.Y,veloData.meanV,npoints); % profile_v is the v velocity in the y direction (measured along x)
    % Subtract bulk flow from v
    veloData.vVelo = veloData.vVelo-mean(veloData.profileVelo);
    fprintf('Mean v flow speed: %f m.s-1\n',mean(veloData.profileVelo))
    elseif opts.flowVeloALreadyKnown == 1
    veloData.vVelo = veloData.vVelo-opts.fixedVelo;
    end


elseif opts.flowDirection == 3
    quiverColor(veloData.X,veloData.Y,veloData.meanU,veloData.meanV,maxScale,vecScale); 
end

end

