%% 
clc
clear
close all

%% GENERATE ELLIPSE 01
imageSizeX = 512;
imageSizeY = 512;
imageSizeZ = 128;
[columnsInImage,rowsInImage,slicesInImage] = meshgrid(1:imageSizeX, 1:imageSizeY, 1:imageSizeZ);
% Next create the ellipse in the image.
centerX = 256;
centerY = 128;
centerZ = 64;

radiusX = 50;
radiusY = 50;
radiusZ = 50;

ellipsePixels = (rowsInImage - centerY).^2 ./ radiusY^2 ...
    + (columnsInImage - centerX).^2 ./ radiusX^2 ...
    + (slicesInImage - centerZ).^2 ./ radiusZ^2 <= 1;


figure(1), 
volshow(ellipsePixels, [1, 1, 1], 0.5, 'r')
title('3D Sphere: Volume 1', 'FontSize', 20);

ro1_rescale = ellipsePixels;
%% GENERATE ELLIPISE 02 
centerX = 256-20;
centerY = 128;
centerZ = 64-10; % Superior

radiusX = 48;
radiusY = 48;
radiusZ = 48;

ellipsePixels = (rowsInImage - centerY).^2 ./ radiusY^2 ...
    + (columnsInImage - centerX).^2 ./ radiusX^2 ...
    + (slicesInImage - centerZ).^2 ./ radiusZ^2 <= 1;


figure(1), hold on
volshow(ellipsePixels, [1, 1, 1], 0.5, 'g'),
title('3D Sphere: Volume 2', 'FontSize', 20);
ro2_rescale = ellipsePixels;
%% VISUALIZATION 

finalSUm = ro1_rescale + ro2_rescale;
figure(10)
for iter1= 1:size(finalSUm , 3)
    figure(10)
    title(iter1)
    tempImg = finalSUm(:,:,iter1);
    imshow(tempImg, [])
    drawnow
end     


%%

localControlDeviationMap100 = calLocalControlDeviation(ro1_rescale,ro2_rescale);

%% TWO-RT_Structure
%
% ro1
targetVolume = double(ro1_rescale);
[r,c,v] = ind2sub(size(targetVolume),find(targetVolume== 1));
orig  = [mean(c) mean(r) mean(v)];

[m,n,p]= size(targetVolume);
[X_temp, Y_temp, Z_temp] =meshgrid(1:n,1:m,1:p);
S = isosurface(X_temp, Y_temp, Z_temp, targetVolume, 0.5 );
faces = S.faces;
vertices = S.vertices;
vert1_1 = vertices(faces(:,1),:);
vert1_2 = vertices(faces(:,2),:);
vert1_3 = vertices(faces(:,3),:);

%rot2
targetVolume2 = double(ro2_rescale);
[m2,n2,p2]= size(targetVolume2);
[X_temp2, Y_temp2, Z_temp2] =meshgrid(1:n2,1:m2,1:p2);
S2 = isosurface(X_temp2, Y_temp2, Z_temp2, targetVolume2, 0.5 );
faces2 = S2.faces;
vertices2 = S2.vertices;
vert2_1 = vertices2(faces2(:,1),:);
vert2_2 = vertices2(faces2(:,2),:);
vert2_3 = vertices2(faces2(:,3),:);


close all
r = 25
figure(21); clf;, hold on, view(3)
volshow(targetVolume, [1, 1, 1], 0.5, 'r'),
volshow(targetVolume2, [1, 1, 1], 0.5, 'g'),
% target volume: ro1 
% trisurf(faces, vertices(:,1),vertices(:,2),vertices(:,3),'FaceAlpha', 0.5)
% trisurf(faces2, vertices2(:,1),vertices2(:,2),vertices2(:,3),'FaceAlpha', 0.9)
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

%%
% Unit of Azimuth, Elevation is radian.
stAzimuth = 0;
endAzimuth = 2*pi;
stElevation = 0;
endElecation = 2*pi;
StepSize = 16;
Stepper = (2*pi)/StepSize;
r = 100;
localControlDeviationMap = []
counterAzimuth = 0;
for iterAzimuth= stAzimuth:Stepper:endAzimuth
    counterAzimuth = counterAzimuth +1
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

        scatter3(dirShift(1), dirShift(2), dirShift(3), 'r*');
        
        % TARGET VOLUME 1
        [intersectTemp,~,~,~,xcoorTemp] = TriangleRayIntersection(orig, sphSet, ...
            vert1_1, vert1_2, vert1_3, 'planeType', 'one sided', 'lineType', 'line','border', 'normal','epsilon', '1e-5'); % 기준이 다름, 명확하게 확인할 것. Direction은 origin에 더해지는 normal vector면 될듯함.
        fprintf('[TARGET VOLUME 1] Number of: faces=%i, points=%i, intresections=%i;\n', ...
            size(faces,1), size(vertices,1), sum(intersectTemp));
%         scatter3(xcoorTemp(intersectTemp,1), xcoorTemp(intersectTemp,2), xcoorTemp(intersectTemp,3), 100, 'b', 'o', 'filled')
        % Line
        %         line('XData',orig(1)+[0 dirShift(1)],'YData',orig(2)+[0 dirShift(2)],'ZData',...
        %   orig(3)+[0 dirShift(3)],'Color','r','LineWidth',3)
        line('XData',[dirShift(1) dirShift_2(1)],'YData',[dirShift(2) dirShift_2(2)],'ZData',...
            [dirShift(3) dirShift_2(3)],'Color','r','LineWidth',3)
        drawnow
        

        % TARGET VOLUME 2
        [intersectTemp2,~,~,~,xcoorTemp2] = TriangleRayIntersection(orig, sphSet, ...
            vert2_1, vert2_2, vert2_3,  'planeType', 'one sided','lineType', 'line', 'border', 'normal', 'epsilon', '1e-5'); % 기준이 다름, 명확하게 확인할 것. Direction은 origin에 더해지는 normal vector면 될듯함.
        fprintf('[TARGET VOLUME 2] Number of: faces=%i, points=%i, intresections=%i;\n', ...
            size(faces2,1), size(vertices2,1), sum(intersectTemp2));
%         scatter3(xcoorTemp2(intersectTemp2,1), xcoorTemp2(intersectTemp2,2), xcoorTemp2(intersectTemp2,3), 200, 'b', 'o', 'filled')
%         % Line
%         %         line('XData',orig(1)+[0 dirShift(1)],'YData',orig(2)+[0 dirShift(2)],'ZData',...
%         %   orig(3)+[0 dirShift(3)],'Color','r','LineWidth',3)
%         line('XData',[dirShift(1) dirShift_2(1)],'YData',[dirShift(2) dirShift_2(2)],'ZData',...
%             [dirShift(3) dirShift_2(3)],'Color','r','LineWidth',3)
        drawnow
        
        % ANALYSIS AND ASSIGN
        
        ro1_points = [xcoorTemp(intersectTemp,1), xcoorTemp(intersectTemp,2), xcoorTemp(intersectTemp,3)];
        ro2_points = [xcoorTemp2(intersectTemp2,1), xcoorTemp2(intersectTemp2,2), xcoorTemp2(intersectTemp2,3)];
        ro1_pointsF = mean(ro1_points, 1);
        ro2_pointsF = mean(ro2_points, 1);
        
%         pointsSet = [ro1_points(1,:);ro2_points(1,:)];  %  
%         pointsSet2 = [ro1_points(2,:);ro2_points(2,:)];  %  
        pointsSetF = [ro1_pointsF; ro2_pointsF];
        
        
        
        if iterElevation > pi
            scatter3(pointsSetF(:,1), pointsSetF(:,2), pointsSetF(:,3), 200, 'b', 'o', 'filled')
            CalculatedDistance1= dist3dSigned(orig, ro1_pointsF, ro2_pointsF);
            localControlDeviationMap(counterElevation, counterAzimuth) = CalculatedDistance1;
        else    
            scatter3(pointsSetF(:,1), pointsSetF(:,2), pointsSetF(:,3), 200, 'k', 'o', 'filled')
            CalculatedDistance2= dist3dSigned(orig,ro1_pointsF, ro2_pointsF);
            localControlDeviationMap(counterElevation, counterAzimuth) = CalculatedDistance2;
        end
        a=0;
    end
end

% localControlDeviationMap(localControlDeviationMap>=50) =NaN;
% localControlDeviationMap(localControlDeviationMap==0) =NaN;
% [localControlDeviationMap_interp,s,exitflag] = smoothn(localControlDeviationMap)

figure, imshow(localControlDeviationMap, [])
ylabel('Elevation')
xlabel('Azimuth')
%%

