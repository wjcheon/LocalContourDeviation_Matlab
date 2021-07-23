function [localControlDeviationMap] = calLocalControlDeviation(reference_,target_)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
visualFlag = false;

% Reference volume
referenceVolume = double(reference_);
[r,c,v] = ind2sub(size(referenceVolume),find(referenceVolume== 1));
orig  = [mean(c) mean(r) mean(v)];

[m,n,p]= size(referenceVolume);
[X_temp, Y_temp, Z_temp] =meshgrid(1:n,1:m,1:p);
S = isosurface(X_temp, Y_temp, Z_temp, referenceVolume, 0.5 );
faces = S.faces;
vertices = S.vertices;
vert1_1 = vertices(faces(:,1),:);
vert1_2 = vertices(faces(:,2),:);
vert1_3 = vertices(faces(:,3),:);

% Target volume
targetVolume = double(target_);
[m2,n2,p2]= size(targetVolume);
[X_temp2, Y_temp2, Z_temp2] =meshgrid(1:n2,1:m2,1:p2);
S2 = isosurface(X_temp2, Y_temp2, Z_temp2, targetVolume, 0.5 );
faces2 = S2.faces;
vertices2 = S2.vertices;
vert2_1 = vertices2(faces2(:,1),:);
vert2_2 = vertices2(faces2(:,2),:);
vert2_3 = vertices2(faces2(:,3),:);

% VISUALIZATION 
if(visualFlag)
    figure(21); clf;, hold on, view(3)
    volshow(referenceVolume, [1, 1, 1], 0.5, 'r'),
    volshow(targetVolume, [1, 1, 1], 0.5, 'g'),
    
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
end

%%
% Unit of Azimuth, Elevation is radian.
stAzimuth = 0;
endAzimuth = 2*pi;
stElevation = 0;
endElecation = 2*pi;
StepSize = 16;
Stepper = (2*pi)/StepSize;
r = 100; % Searching radius for intercept points 
localControlDeviationMap = [];
counterAzimuth = 0;
f = waitbar(0,'Please wait...');
for iterAzimuth= stAzimuth:Stepper:endAzimuth
    counterAzimuth = counterAzimuth +1;
    counterElevation = 0;
    for iterElevation= stElevation:Stepper:endElecation
        counterElevation = counterElevation+1;
        % LINE DEFINITION
        azimuthTemp = iterAzimuth;
        elevationTemp = iterElevation;
        [xTemp, yTemp, zTemp] = sph2cart(azimuthTemp,elevationTemp,r);
        
        sphSet(1) = xTemp;
        sphSet(2) = yTemp;
        sphSet(3) = zTemp;
        
        dirShift(1) = orig(1)+xTemp;
        dirShift(2) = orig(2)+yTemp;
        dirShift(3) = orig(3)+zTemp;
        
        dirShift_2(1) = orig(1)-xTemp;
        dirShift_2(2) = orig(2)-yTemp;
        dirShift_2(3) = orig(3)-zTemp;

%         scatter3(dirShift(1), dirShift(2), dirShift(3), 'r*');
        
        % Calculate an intercept point of reference volume 
        [intersectTemp,~,~,~,xcoorTemp] = TriangleRayIntersection(orig, sphSet, ...
            vert1_1, vert1_2, vert1_3, 'planeType', 'one sided', 'lineType', 'line','border', 'normal','epsilon', '1e-5'); % 기준이 다름, 명확하게 확인할 것. Direction은 origin에 더해지는 normal vector면 될듯함.
%         fprintf('[TARGET VOLUME 1] Number of: faces=%i, points=%i, intresections=%i;\n', ...
%             size(faces,1), size(vertices,1), sum(intersectTemp));

        

        % Calculate an intercept point of target volume
        [intersectTemp2,~,~,~,xcoorTemp2] = TriangleRayIntersection(orig, sphSet, ...
            vert2_1, vert2_2, vert2_3,  'planeType', 'one sided','lineType', 'line', 'border', 'normal', 'epsilon', '1e-5'); % 기준이 다름, 명확하게 확인할 것. Direction은 origin에 더해지는 normal vector면 될듯함.
%         fprintf('[TARGET VOLUME 2] Number of: faces=%i, points=%i, intresections=%i;\n', ...
%             size(faces2,1), size(vertices2,1), sum(intersectTemp2));

        
        % ANALYSIS AND ASSIGN
        ro1_points = [xcoorTemp(intersectTemp,1), xcoorTemp(intersectTemp,2), xcoorTemp(intersectTemp,3)];
        ro2_points = [xcoorTemp2(intersectTemp2,1), xcoorTemp2(intersectTemp2,2), xcoorTemp2(intersectTemp2,3)];
        ro1_pointsF = mean(ro1_points, 1);
        ro2_pointsF = mean(ro2_points, 1);
        pointsSetF = [ro1_pointsF; ro2_pointsF];
        
        
        if iterElevation > pi
            CalculatedDistance1= dist3dSigned(orig, ro1_pointsF, ro2_pointsF);
            localControlDeviationMap(counterElevation, counterAzimuth) = CalculatedDistance1;
            if(visualFlag)
               scatter3(pointsSetF(:,1), pointsSetF(:,2), pointsSetF(:,3), 200, 'b', 'o', 'filled') 
            end
            
        else    
            CalculatedDistance2= dist3dSigned(orig,ro1_pointsF, ro2_pointsF);
            localControlDeviationMap(counterElevation, counterAzimuth) = CalculatedDistance2;
            if(visualFlag)
               scatter3(pointsSetF(:,1), pointsSetF(:,2), pointsSetF(:,3), 200, 'k', 'o', 'filled') 
            end
        end
    end
    waitbar(iterAzimuth/endAzimuth,f,'Processing your data');
end
waitbar(iterAzimuth/endAzimuth,f,'Finishing');
close(f)

% Smoothing for 'NaN' value.
% localControlDeviationMap(localControlDeviationMap>=50) =NaN;
% localControlDeviationMap(localControlDeviationMap==0) =NaN;
% When there is no "an intercept point", NaN value was returend.
[localControlDeviationMap,s,exitflag] = smoothn(localControlDeviationMap)

figure, imshow(localControlDeviationMap, [])
ylabel('Elevation')
xlabel('Azimuth')

end

