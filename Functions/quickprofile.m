function [profile_u,profile_r,profile_x,profile_y,xp,yp] = quickprofile(X,Y,velo,npoints) % meshgrid X, Y. Velo is the corresponding u or v data along the path 
% Create a quick velocity profile along a line in a velocity field.
% If points are not defined, get them from a plot using ginput.

% Options
plotOn = 1;

% query two points along the horizontal or vertical axes to measure the
% bulk flow
[xp,yp] = ginput(2);
% hold on
plot(xp(1),yp(1),'go')
plot(xp(2),yp(2),'ro')
plot(xp,yp,'k--')


profile_x = linspace(xp(1),xp(2),npoints)'; % x coordinates of the interpolated points along the line
profile_y = linspace(yp(1),yp(2),npoints)'; % y coordinates of the interpolated points along the line
profile_r = sqrt((profile_x - profile_x(1)).^2 +...
                 (profile_y - profile_y(1)).^2); % enables plotting alng one line
profile_u = interp2(X,Y,velo,profile_x,profile_y); % interpolate along the line


if plotOn == 1
    figure
    plot(profile_r,profile_u,'.')
    xlabel('length of profile (m)')
    ylabel('velocity (m s^{-1})')

end
