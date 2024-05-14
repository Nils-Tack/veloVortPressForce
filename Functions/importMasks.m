function BW = importMasks(D_BW,f)
if nargin > 1
% Import image data for one mask # f
BW = importdata(quickfilepath(D_BW(f)));

else % import the full stack of BW and store in a 3D array
    for i = 1:length(D_BW)
        BW(:,:,i) = importdata(quickfilepath(D_BW(i)));
    end

end
end