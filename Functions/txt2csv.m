function txt2csv(D,path)
% Options
filestub = 'B';         % Default file name in Queen2
opts.convertUnits = 1;  % Option to enable (1) or disable (0) unit conversion

fprintf('PIV files conversion...');
 for i=1:length(D.PIVraw)
     % Generate name of PIVclean files
     num = num2str(i,'%05d'); % create file number with 5 digits
     rf = strcat(filestub, num); % concatenate file number; filestub + file number
     
     % Import raw data
     velo = readmatrix(quickfilepath(D.PIVraw(i))); % read data

     % Convert units if necessary
     if opts.convertUnits == 1
     velo(:,1:2) = velo(:,1:2)/1000; % !!!! convert mm to m because DaVis exports x and y data in mm!!!!
     end

     % Export clean and converted data
     writematrix(velo,fullfile(path.PIVclean,[sprintf('B%05d',i),'.csv'])); % export clean data; accepted format .csv

     progressCount(i,length(D.PIVraw)); % display export progress
end

fprintf('done\n');

end