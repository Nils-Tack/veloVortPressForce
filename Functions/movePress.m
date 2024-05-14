function movePress(D,path)
fprintf('Moving press files...');

for id = 1:length(D.preliminaryPress)
    filepath_in = fullfile(path.PIVclean,D.preliminaryPress(id).name);
    filepath_out = fullfile(path.pressure,D.preliminaryPress(id).name);
    movefile(filepath_in, filepath_out);
end

fprintf('done\n');

end