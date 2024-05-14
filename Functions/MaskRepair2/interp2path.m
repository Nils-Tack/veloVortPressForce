function [pts_out,s_in,s_out] = interp2path(pts_in,interp_method,norm_n)
% Interpolate 2D data along a path
% coords_in is an array with two columns (x,y) ordered by distance along a path
% n is number of points to return
% interp_method is interpolation method to use, e.g., linear, pchip
% norm_n: 0 - return n points, 1 - return n times path length points

% Get distance along path
x_in = pts_in(:,1);
y_in = pts_in(:,2);
s_in = zeros(size(x_in,1),1);
s_in(2:end) = cumsum(sqrt( diff(x_in).^2 + diff(y_in).^2 ));

% Normalize number of points to interpolate over to path length
if norm_n
    distBetweenPt = 0.0005; % desired distance between interpolated points; in meter
    n = floor(s_in(end)/distBetweenPt);
   % n = floor(n*s_in(end));
end

% Interpolate
s_out = linspace(0,s_in(end),n)';
x_out = interp1(s_in,x_in,s_out,interp_method);
y_out = interp1(s_in,y_in,s_out,interp_method);
pts_out = [x_out y_out];

end