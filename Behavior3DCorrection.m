function Behavior3DCorrection
%{Behavior3DCorrection is designed to merge and correct x,z (side view) and
%x,y (top view) tracking data derived from video of an organism (a squid 
%paralarvae) swimming in a cubic arena taken with two perpendicular cameras.
%This code assumes the data is stored in a single Excel spreadsheet where 
%the columns are: time (s), xs, zs, xt, yt. These known data are input into 
%an overconstrained system of equations and assessed using least sum of 
%squares to approximate the true x, y, and z values with the most minimal 
%error. This code accompanies the publication: Zakroff, C. Mooney, TA, 
%Wirth, C. Ocean Acidification Responses in Paralarval Squid Swimming
%Behavior Using a Novel 3D Tracking System. (2017). Hydrobiologia Vol: Pages.
%DOI:
%
%Version 1.5 written by Casey Zakroff (czakroff@whoi.edu) May 11 2017
%in MATLAB version 2016b on Mac. Code and protocols available at: 
%https://github.com/czakroff/3D-Swimming-Behavior
%}

%User Input of Excel Data Sheet
filename = input('\nEnter the name of the Excel file: ', 's');
sheet = input('Enter the name of the Excel sheet: ', 's');
dataRange = input('Enter the range of the data cells, i.e. "A5:E3602":', 's'); %This data range works for 2 minute videos (3598 frames)
pLcode = strcat(filename,'_',sheet); %ID code for each individual squid paralarvae (pL)

%Read in t, xs, zs, xt, yt
rawData = xlsread(filename, sheet, dataRange); 

%Establish known values of from your arena
Q = 9.6; %Q is the length of the side of your cubic arena (in cm)
%Read in Qbs - the length of the back of the arena in the side video image
Qbs = input('Enter Qbs value: ');
%Read in Qbt - the length of the back of the arena in the top video image
Qbt = input('Enter Qbt value: '); 

%Calculate the hypoteneuse of the right triangle between Q and Q back for both the side and top.
qbs = (Q-Qbs)/2*sqrt(2);
qbt = (Q-Qbt)/2*sqrt(2);

%Set up output data array and assign time to first column
correctedData = zeros(size(rawData,1),4);
correctedData(:,1) = rawData(:,1);

%Assign xs, zs, xt, and yt
xs = rawData(1:size(rawData,1),2);
zs = rawData(1:size(rawData,1),3);
xt = rawData(1:size(rawData,1),4);
yt = rawData(1:size(rawData,1),5);

%Assign limits and options for fmincon
A = ones(1,11)*0.001;
B = ones(1,11)*Q;
options = optimoptions('fmincon','Display','final');

for i = 1:size(rawData,1)
    %Assign guesses for fmincon
    % x0 = [x, y, z, xis, zis, xit, yit, qis, qit, Qs, Qt]
    x0 = [xs(i), yt(i), zs(i), xs(i), zs(i), xt(i), yt(i), 0, 0, Q, Q];
    
    %Run least sum of squares using fmincon
    [p,fval,exitflag,output] = fmincon(@leastSumSqr,x0,[],[],[],[],A,B,[],options);
    correctedData(i,2:4) = p(1:3); %store corrected t,x,y,z data
end

%Output corrected data. Creates CSV on Macs
xlswrite(pLcode, correctedData);

%Call analysis code to get 3D metrics and basic plots
Behavior3DAnalysis(correctedData,pLcode);

%System of equations assessed by fmincon
    function fVal = leastSumSqr(x0)
    
        F = zeros(1,12);
        
        F(1) = x0(1)*x0(10)/x0(4)-Q;
        F(2) = x0(3)*x0(10)/x0(5)-Q;
        F(3) = x0(1)*x0(11)/x0(6)-Q;
        F(4) = x0(2)*x0(11)/x0(7)-Q;
        F(5) = (x0(4)+x0(8)/sqrt(2))-xs(i);
        F(6) = (x0(5)+x0(8)/sqrt(2))-zs(i);
        F(7) = (x0(6)+x0(9)/sqrt(2))-xt(i);
        F(8) = (x0(7)+x0(9)/sqrt(2))-yt(i);
        F(9) = (x0(10)+2*x0(8)/sqrt(2))-Q;
        F(10) = (x0(11)+2*x0(9)/sqrt(2))-Q;
        F(11) = x0(8)/x0(2)-qbs/Q;
        F(12) = x0(9)/x0(3)-qbt/Q;
        
        fVal = sum(F.^2);
    end

end

