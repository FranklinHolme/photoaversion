function [ deltaA ] = deltaAngle(a, b)
% Computes smallest difference between two angles based on this code
% snippet from stack overflow: 
% a = targetA - sourceA
% a -= 360 if a > 180
% a += 360 if a < -180

    deltaA = a - b;
    
    if deltaA > 180
        deltaA = deltaA - 360;
    end

    if deltaA < -180
        deltaA = deltaA + 360;
    end


end

