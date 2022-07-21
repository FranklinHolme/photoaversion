% Compute speed from arrays of x and y positions. Returns speed in
% distance/s
function [speed] = positionToSpeed(X, Y, framerate)
    speed = NaN * ones(size(X));
    for i = 2:length(X)
        speed(i) = sqrt( (X(i) - X(i-1))^2 + (Y(i) - Y(i-1))^2);
    end 
    speed = speed * framerate; 
end