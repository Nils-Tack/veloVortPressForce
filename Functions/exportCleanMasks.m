function exportCleanMasks(path,D,sizeBW,opts)

filepath_masks = fullfile(path.masksBW,'Clean masks');
mkdir(filepath_masks)

fprintf('Exporting BW masks from outlines...')
for i = 1:length(D.outlines)
outline = importdata(quickfilepath(D.outlines(i)));
BWoutline = flipud(poly2mask(round(outline(:,1)*opts.scale),round(outline(:,2)*opts.scale),sizeBW(1),sizeBW(2)));

filenameBW = sprintf('BW_%05g',i);
imwrite(BWoutline,fullfile(filepath_masks,[filenameBW, '.tif']))
progressCount(i,length(D.outlines)); % display export progress
end

fprintf('done\n');
end