function AllData = importData(D_data)

 for i = 1:length(D_data)
        AllData(:,:,i) = importdata(quickfilepath(D_data(i)));
 end

end