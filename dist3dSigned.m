
function [dist] = dist3dSigned(origin,P1,P2)

    distanceTemp1 = dist3d(P1,origin);
    distanceTemp2 = dist3d(P2,origin);
    dist = distanceTemp1-distanceTemp2; % Reference - Target 
 
end