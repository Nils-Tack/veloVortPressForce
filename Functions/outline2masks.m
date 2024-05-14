function BW = outline2masks(D,opts)
imageSize = size(importMasks(D.masksBW,1));    %Measure the length and width in pixel of one original BW mask.

fprintf('Converting outlines to BW masks...');
for i = 1:length(D.outlines)
    outline = importdata(quickfilepath(D.outlines(i)))*opts.scale; % load the corresponding outline and convert it to pixel based on the scale
    BW(:,:,i) = poly2mask(outline(:,1),outline(:,2),imageSize(1),imageSize(2)); % transform the outline to BW mask, using the same length and width as the original mask.
    progressCount(i,length(D.outlines)); % display export progress
end
fprintf('done\n');

end