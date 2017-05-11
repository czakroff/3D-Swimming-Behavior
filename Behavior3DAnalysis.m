function Behavior3DAnalysis(data, filename)
%{Behavior3DAnalysis takes the corrected x, y, z positional data from 
%Behavior3DCorrection and calculates total distance traveled, average,
%peak, min, max, and median of 3D, horizontal, and vertical velocity, 3D
%volume covered, 2D turning angles, and path tortuosity (a ratio of the 
%distance transitited between two points and the true distance between those
%points) of an individual swimming organism. These data are output to a 
%'results' text file and to an Excel file named with the individual 
%organism ID code. A 3D plot of the organism's swimming track, a plot of 
%the tortuosity of that path, and a 3D polygon of the volume covered by the
%organism are produced and output as tiff files. A histogram of the
%organism's turning angles over the entire track is output as a png.
%This code accompanies the publication: Zakroff, C. Mooney, TA, 
%Wirth, C. Ocean Acidification Responses in Paralarval Squid Swimming
%Behavior Using a Novel 3D Tracking System. (2017). Hydrobiologia. Vol: Pages.
%DOI:
%
%Version 1.5 written by Casey Zakroff (czakroff@whoi.edu) May 11 2017
%in MATLAB version 2016b on Mac. Code and protocols available at: 
%https://github.com/czakroff/3D-Swimming-Behavior
%}

%Set variables to set sampling resolution for tortuosity and turning angle
res = 30; %Sampling resolution (value is in frames; 30 fps = 1 s segments)
numPts = floor(size(data,1)/res); %Number of points sampled from the data

%Create arrays to store 3D metrics
distance = zeros(size(data,1)-1,1); 
vel = zeros(size(data,1)-1,1); %3D velocity
vertVel = zeros(size(data,1)-1,1); %vertical velocity
horVel = zeros(size(data,1)-1,1); %horizontal velocity
angles = zeros(numPts-1,1); %turning angles
tort = zeros(size(data,1)-res-1,1); %tortuosity

%Calculate total distance traveled (cm) and velocities (cm/s)
%by looping through all data and recording values between successive 
%positions.
for i = 1:(size(data,1)-1)
    tDiff = data(i+1,1)-data((i),1); %difference in time (s)
    xDiff = data(i+1,2)-data((i),2); %difference in x (cm)
    yDiff = data(i+1,3)-data((i),3); %difference in y (cm)
    zDiff = data(i+1,4)-data((i),4); %difference in z (cm)
    distance(i) = sqrt(xDiff^2+yDiff^2+zDiff^2);
    vel(i) = distance(i)/tDiff;
    vertVel(i) = zDiff/tDiff;
    horVel(i) = sqrt(xDiff^2+yDiff^2)/tDiff;
end

%Store distance and speed data
totDist = sum(distance);
avgVel = mean(vel);
peakVel = max(vel);
minVel = min(vel);
medVel = median(vel);
avgVVel = mean(vertVel);
peakVVel = max(vertVel);
minVVel = min(vertVel);
medVVel = median(vertVel);
avgHVel = mean(horVel);
peakHVel = max(horVel);
minHVel = min(horVel);
medHVel = median(horVel);

%Calculate hull (for polygon) and 3D volume covered
[K,vol3D] = convhull(data(:,2),data(:,3),data(:,4));

%Calculate 2D turning angles between successive 1-second steps (res)
j = 1;
for i = 1:(numPts-1)
    xdiff1 = data((j+res),2)-data(j,2);
    ydiff1 = data((j+res),3)-data(j,3);
    xdiff2 = data((j+2*res),2)-data((j+res),2);
    ydiff2 = data((j+2*res),3)-data((j+res),3);
    a = [xdiff1,ydiff1];
    b = [xdiff2,ydiff2];
    angle = atan2(abs(det([a;b])),dot(a,b));
    angles(i) = angle*180/pi;   
    j = j+res;
end

%Store turning angle data
avgAngl = mean(angles);
maxAngl = max(angles);
minAngl = min(angles);
medAngl = median(angles);

%Tortuosity
%Loop for calculating tortuosity over 1-second steps (res)
for i = 1:(size(data,1)-res-1)
    xDiffT = data((i+res),2)-data((i),2);
    yDiffT = data((i+res),3)-data((i),3);
    zDiffT = data((i+res),4)-data((i),4);
    distPts = sqrt(xDiffT^2+yDiffT^2+zDiffT^2); %distance between points
    distTort = sum(distance(i:(i+res))); %distance covered by organism
    tort(i) = distTort/distPts;
end

%Store tortuosity data
avgTort = mean(tort);
peakTort = max(tort);
minTort = min(tort);
medTort = median(tort);

figure
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPositionMode = 'manual';
fig.PaperSize = [33.8 33.8];

%Plot 3D track and output tiff
plot3(data(:,2),data(:,3),data(:,4)-9.6);
ax = gca;
ax.FontUnits = 'centimeters';
ax.FontName = 'Times New Roman';
ax.FontSize = 0.35;
ax.XTickMode = 'manual';
ax.YTickMode = 'manual';
ax.XLimMode = 'manual';
ax.YLimMode = 'manual';
ax.XLim = [0,9.600];
ax.YLim = [0,9.600];
ax.ZLim = [-9.600,0];
ax.XTick = [0,1,2,3,4,5,6,7,8,9];
ax.XTickLabel = {' ', '1', ' ', '3', ' ', '5', ' ', '7', ' ', '9'};
ax.YTick = [0,1,2,3,4,5,6,7,8,9];
ax.YTickLabel = {' ', '1', ' ', '3', ' ', '5', ' ', '7', ' ', '9'};
ax.ZTick = [-9,-8,-7,-6,-5,-4,-3,-2,-1,0];
ax.ZTickLabel = {'-9', '', '-7', '', '-5', '', '-3', '', '-1', ''};
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.ZGrid = 'on';
ax.XLabel.String = 'Length (cm)';
ax.YLabel.String = 'Length (cm)';
ax.ZLabel.String = 'Depth (cm)';
print('-dtiff','-r600',filename);

%Calculate colormap for tortuosity plot
[t1,t2] = size(tort);
[d1,d2] = size(data(:,2));
col = tort;
for i = t1:d1-1
    col = [col; tort(t1)];
end

%Plot tortuosity as thin 3D surface and output tiff
h = surface([data(:,2), data(:,2)], [data(:,3), data(:,3)],...
    [data(:,4)-9.6, data(:,4)-9.6],[col,col],...
    'FaceColor','none','EdgeColor','interp','LineWidth',2);
colormap(jet);
caxis([1 5]);
c = colorbar;
text('String','Tortuosity', 'Rotation', 90,...
    'Position', [242.25 302.5 -227.2],...
    'FontUnits', 'centimeters', 'FontSize', 0.35,...
    'FontName', 'Times New Roman', 'color','k');
print('-dtiff','-r400',strcat(filename,'_tortuosity'));

%Create 3D volume polygon and output tiff
colormap(parula);
trisurf(K,data(:,2),data(:,3),data(:,4)-9.6);
ax = gca;
ax.FontUnits = 'centimeters';
ax.FontName = 'Times New Roman';
ax.FontSize = 0.35;
ax.XTickMode = 'manual';
ax.YTickMode = 'manual';
ax.XLimMode = 'manual';
ax.YLimMode = 'manual';
ax.XLim = [0,9.600];
ax.YLim = [0,9.600];
ax.ZLim = [-9.600,0];
ax.XTick = [0,1,2,3,4,5,6,7,8,9];
ax.XTickLabel = {' ', '1', ' ', '3', ' ', '5', ' ', '7', ' ', '9'};
ax.YTick = [0,1,2,3,4,5,6,7,8,9];
ax.YTickLabel = {' ', '1', ' ', '3', ' ', '5', ' ', '7', ' ', '9'};
ax.ZTick = [-9,-8,-7,-6,-5,-4,-3,-2,-1,0];
ax.ZTickLabel = {'-9', '', '-7', '', '-5', '', '-3', '', '-1', ''};
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.ZGrid = 'on';
ax.XLabel.String = 'Length (cm)';
ax.YLabel.String = 'Length (cm)';
ax.ZLabel.String = 'Depth (cm)';
print('-dtiff','-r600',strcat(filename,'_volume'));

%Plot histogram of turning angles and output png
angBins = 18; %Number of bins for angle histogram
histogram(angles,angBins,'BinWidth', 10,'Normalization','probability');
ax = gca;
ax.XTickMode = 'manual';
ax.XLimMode = 'manual';
ax.FontUnits = 'centimeters';
ax.FontName = 'Times New Roman';
ax.FontSize = 0.35;
ax.Box = 'Off';
ax.XLabel.String = 'Turning Angle (degrees)';
ax.XLim = [0 180];
ax.YLim = [0 0.5];
ax.XTick = [0,20,40,60,80,100,120,140,160,180];
ax.YTick = [0,0.1,0.2,0.3,0.4,0.5];
ax.YLabel.String = 'Frequency';
print(strcat(filename,'_turningAngles'), '-dpng');
    
%Compile results
results = [totDist; avgVel; peakVel; minVel; medVel; avgVVel; 
    peakVVel; minVVel; medVVel; avgHVel; peakHVel; minHVel; medHVel;
    vol3D; avgAngl; maxAngl; minAngl; medAngl; avgTort; peakTort; minTort
    medTort];

%Output data to results textfile
fid = fopen('results.txt', 'a');
spec = '\n\n%s%123.4f%123.4f%123.4f%123.4f%123.4f%123.4f%123.4f%123.4f%123.4f%123.4f%123.4f%123.4f%123.4f%123.4f%123.4f%123.4f%123.4f%123.4f%123.4f%123.4f%123.4f%123.4f';
fprintf(fid, spec, filename, results);
fclose(fid);

%{
%Output data on PC (output as one Excel file)
xlswrite(filename,results,'results');
xlswrite(filename,distance, 'distance');
xlswrite(filename,vel,'3Dvelocity');
xlswrite(filename,vertVel,'verticalVelocity');
xlswrite(filename,horVel,'horizontalVelocity');
xlswrite(filename,angles, 'angles');
xlswrite(filename,tort, 'tortuosity');
%}

%Output data on Mac (output as individual CSVs)
xlswrite(strcat(filename,'_results'),results);
xlswrite(strcat(filename,'_distance'),distance);
xlswrite(strcat(filename,'_3Dvelocity'),vel);
xlswrite(strcat(filename,'_verticalVelocity'),vertVel);
xlswrite(strcat(filename,'_horizontalVelocity'),horVel);
xlswrite(strcat(filename,'_angles'),angles);
xlswrite(strcat(filename,'_tortuosity'),tort);

end

