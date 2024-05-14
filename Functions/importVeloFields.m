function veloData = importVeloFields(D)

% Initiate X, Y, u and v (pre-allocate the variable for speed); store in
% one structure
tempVelo = importdata(quickfilepath(D.PIVclean(1)));
veloData.X = reshape(tempVelo(:,1),length(unique(tempVelo(:,1))),length(tempVelo(:,1))/length(unique(tempVelo(:,1))))';
veloData.Y = reshape(tempVelo(:,2),length(unique(tempVelo(:,1))),length(tempVelo(:,1))/length(unique(tempVelo(:,1))))';
veloData.uVelo = zeros(size(veloData.X,1),size(veloData.X,2),length(D.PIVclean));
veloData.vVelo = zeros(size(veloData.X,1),size(veloData.X,2),length(D.PIVclean));
veloData.UVspeed = zeros(size(veloData.X,1),size(veloData.X,2),length(D.PIVclean));

fprintf('Velocity files extraction...');
    for i = 1:length(D.PIVclean)
        tempVeloData = importdata(quickfilepath(D.PIVclean(i)));
        veloData.uVelo(:,:,i) = reshape(tempVeloData(:,3),length(unique(tempVeloData(:,1))),length(tempVeloData(:,1))/length(unique(tempVeloData(:,1))))';
        veloData.vVelo(:,:,i) = reshape(tempVeloData(:,4),length(unique(tempVeloData(:,1))),length(tempVeloData(:,1))/length(unique(tempVeloData(:,1))))';
        veloData.UVspeed(:,:,i) = sqrt(veloData.uVelo(:,:,i).^2 + veloData.vVelo(:,:,i).^2); % compute flow velocity magnitude
        progressCount(i,length(D.PIVclean)); % display export progress
    end
fprintf('done\n');

end