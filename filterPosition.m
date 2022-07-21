% Filters position data consisting of x,y coordinate pairs with a small
% median filter and a big mean filter. 
function [position_out] = filterPosition(position_in, framerate)
    position_out = medfilt1(position_in, round(framerate / 5), 'omitnan'); % Third argument leaves out the NaN values when calculating median 
    position_out = movmean(position_out, round(framerate / 2), 'omitnan'); 
end