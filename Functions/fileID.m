function filename = fileID(D_files, f)
if nargin > 1
% display title of the figure as file #
temp = split(D_files(f).name,'.'); % change delimiter accordingly
filename = temp{1};

else
% store the name of each file in a stack to use as individual titles in figures.
    for i = 1:length(D_files)
    temp = split(D_files(i).name,'.'); % change delimiter accordingly
    filename{i,1} = temp{1}; % store names in array
    end
end


end