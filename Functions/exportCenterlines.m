function exportCenterlines(path,Data)
% Make centerline folder in te main Data folder
filepathOut = fullfile(path.main,'centerlines');
filepathHeading = fullfile(filepathOut,'heading');
mkdir(filepathOut)
mkdir(filepathHeading)

% export individual centerlines to a .csv file
for i = 1:size(Data.centerlinesMeter,3)
filenameOutline = sprintf('centerline_%05g', i); % Number sequentially
writematrix(Data.centerlinesMeter(:,:,i),fullfile(filepathOut,[filenameOutline,'.csv']))
progressCount(i,size(Data.centerlinesMeter,3)); % display export progress
end

% export the fish headings (only one file)
writematrix(Data.heading,fullfile(filepathHeading,['swimming direction','.csv']))

end