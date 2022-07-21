function [a_filt] = filterAngleData(a, sampling_rate)
% Filters angle data, first with a median filter and then with a mean
% filter. The function uses sampling_rate (in Hz) to set the order (how
% many samples to use before and after a single sample to set its new
% value). This code implements the function from this wikipedia page: 
% https://en.wikipedia.org/wiki/Mean_of_circular_quantities
    
    % Initialize the output array 
    a_filt = NaN * ones(size(a));
    
    % This filtering does not work on angles that are exactly zero or 360. 
    % Add a tiny amount to every angle to make sure these values don't
    % occur
    a = a + 0.0001;

    % Take the sine and cosine of the angle vector (provided in degrees) 
    cosa = cosd(a);
    sina = sind(a);
    
    % Apply a median filter. Second argument is filter order
    cosa_filt = medfilt1(cosa, round(sampling_rate/5), 'omitnan'); 
    sina_filt = medfilt1(sina, round(sampling_rate/5), 'omitnan');
    
    % Apply a mean filter 
    cosa_filt = movmean(cosa_filt, round(sampling_rate/2), 'omitnan');
    sina_filt = movmean(sina_filt, round(sampling_rate/2), 'omitnan'); 
    
    % Reconstruct the filtered sine and cosines of the angle data into
    % angle data 
    arc_tan = atand(sina_filt ./ cosa_filt);
    a_filt(cosa_filt > 0 & sina_filt > 0) = arc_tan(cosa_filt > 0 & sina_filt > 0);
    a_filt(cosa_filt < 0) = arc_tan(cosa_filt < 0) + 180;
    a_filt(cosa_filt > 0 & sina_filt < 0) = arc_tan(cosa_filt > 0 & sina_filt < 0) + 360; 

end