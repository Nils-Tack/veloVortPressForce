function pressData = importPressFields(D)

% Initiate X, Y, and P (pre-allocate the variable for speed); store in
% one structure
tempPress = importdata(quickfilepath(D.pressure(1)));
[pressData.X,pressData.Y] = meshgrid(unique(tempPress(:,1),'stable'),unique(tempPress(:,2),'stable'));
pressData.P = zeros(size(pressData.X,1),size(pressData.X,2),length(D.pressure));

fprintf('Pressure files extraction...');
    for i = 1:length(D.pressure)
        tempPressData = importdata(quickfilepath(D.pressure(i)));
        pressData.P(:,:,i) = reshape(tempPressData(:,7),[size(pressData.X,1),size(pressData.X,2)]);
        % pressData.P(:,:,i) = reshape(tempPressData(:,7),length(unique(tempPressData(:,1))),length(tempPressData(:,1))/length(unique(tempPressData(:,1))))';
        progressCount(i,length(D.pressure)); % display export progress
    end

fprintf('done\n');

end