function [g, NR, SI, TI] = regiongrow(f, S, T)
% REGIONGROW perfom segmentation by region growing.
% [G, NR, SI, TI] = REGIONGROW(F, SR, T), S can be an array (the same size
% as F) with a 1 at the coordinates of every seed point and 0s elsewhere. S
% can also be a single seed value. Similary, T can be and array (the same
% size as F) containing a threshold value for each pixel in F. T can also
% be a scalar, in which case it becomes a global threshold. On the output,
% G is the result of region growing, with each region labeled by a
% different integer, NR is the number of regions, SI is the final seed
% image used by the algorithm, and TI is the image consisting of the pixels
% in F that satisfied the threshold test. 

f = double(f);
% If S is a scalar, obtain the seed image. 
if numel(S) == 1
    SI = f == S;
    S1 = S;
else 
    % S is an array. Eliminate duplicate, connected seed locations to
    % reduce the number of loop executions in the following sections of
    % code.
    SI = bwmorph(S,'shrink',Inf);
    J = find(SI);
    S1 = f(J); % Array of seed values.
end

TI = false(size(f));
for K =1:length(S1)
    seedvalue = S1(K);
    S = abs(f-seedvalue) <= T;
    TI = TI | S;
end
% Use function imreconstruct with SI as the marker image to obtain the
% regions corresponding to each seed in S. Function bwlabel assings a
% different integer to each connected region.
[g, NR] = bwlabel(imreconstruct(SI,TI));
end