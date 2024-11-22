clearvars

%% Defining Directory Location

csvDirList = {'Analysis/Kir2.1 No Escape';
    'Analysis/Kir2.1 Escape';
    'Analysis/DL(Wild) No Escape';
    'Analysis/DL(Wild) Escape'};
csvDirList2 = {'Box Analysis/Kir2.1 No Escape';
    'Box Analysis/Kir2.1 Escape';
    'Box Analysis/DL(Wild) No Escape';
    'Box Analysis/DL(Wild) Escape'};
xlsDir = 'Analysis/Time_fps';
xlsDir2 = 'Analysis/Time of Contact';

%%

maxSpeed = cell(1,4);
avgSpeed = cell(1,4);
peakAcc = cell(1,4);
accAvg = cell(1,4);
accAvg2 = cell(1,4);
accMed2 = cell(1,4);
sizeSpeed = cell(1,4);
initialSizeSpeed = cell(1,4);
avgSizeSpeed = cell(1,4);
angSize = cell(1,4);
timeOfCont = cell(1,4);
distTOC = cell(1,4);
tvVisualAng = cell(1,4);
avgAzimuth = cell(1,4);
avgElevation = cell(1,4);
Type = struct();

for p = 1:4 % look through each of the four cases: Kir No Escape, Kir Escape, DL No Escape, DL Escape
    
    csvDir = fullfile(pwd,csvDirList{p}); %Tracking information
    csvFolders = dir(csvDir);
    csvFolders = {csvFolders(3:end).name}; %Removes invisible folders at beginning of directory and only takes the names
    
    csvDir2 = csvDirList2{p}; %Fly axis information (for orientation)
    csvFolders2 = dir(csvDir2);
    csvFolders2 = {csvFolders2(3:end).name};
    
    xlsFolders = dir(xlsDir); %Time matrix information
    xlsFolders = {xlsFolders(3).name};
    
    xlsFolders2 = dir(xlsDir2); %Time of contact/escape information
    xlsFolders2 = {xlsFolders2(3).name};
    
    maxSpeed{p} = zeros(numel(csvFolders),1);
    avgSpeed{p} = zeros(numel(csvFolders),1);
    accAvg{p} = zeros(numel(csvFolders),1);
    accAvg2{p} = zeros(numel(csvFolders),1);
    accMed2{p} = zeros(numel(csvFolders),1);
    sizeSpeed{p} = zeros(numel(csvFolders),1);
    initialSizeSpeed{p} = zeros(numel(csvFolders),1);
    avgSizeSpeed{p} = zeros(numel(csvFolders),1);
    angSize{p} = zeros(numel(csvFolders),1);
    timeOfCont{p} = zeros(numel(csvFolders),1);
    distTOC{p} = zeros(numel(csvFolders),1);
    tvVisualAng{p} = zeros(numel(csvFolders),1);
    peakAcc{p} = zeros(numel(csvFolders),1);
    avgAzimuth{p} = zeros(numel(csvFolders),1);
    avgElevation{p} = zeros(numel(csvFolders),1);
    
    
    
    for d = 1:numel(csvFolders) % Loops through each predation attempt in each case
        %% STEP 1: Bringing in X, Y, Z values of Damselfly Head/Tail and Fruitfly (Including Axis Coordinates)
        
        csvPath = fullfile(csvDir,csvFolders{d},[csvFolders{d} '_data_xyzpts.csv']); %path for XYZ points of damselfly and fly
        csvPath2 = fullfile(csvDir2,csvFolders2{d},[csvFolders2{d} '_data_xyzpts.csv']); %path for points for fly orientation
        
        Type(p).Data(d).path = csvPath;
        
        XYZ = csvread(csvPath,1,0); %Read in XYZ points from csv
        
        Type(p).Data(d).damXYZ = XYZ(:,1:3); %Damselfly head XYZ data
        damXYZ = Type(p).Data(d).damXYZ;
        
        Type(p).Data(d).flyXYZ = XYZ(:,7:9); %Fly centroid XYZ data
        flyXYZ = Type(p).Data(d).flyXYZ;
        
        Type(p).Data(d).damtail = XYZ(:,4:6); %Damselfly tail XYZ data, does not get used in remainder of script
        
        FlyAxis = csvread(csvPath2,1,0); %Read in fly axis points for calculating fly orientation vector
        Type(p).Data(d).head = FlyAxis(1,1:3); %fly head
        head = Type(p).Data(d).head;
        
        Type(p).Data(d).tail = FlyAxis(1,4:6);
        tail = Type(p).Data(d).tail; %fly tail
        
        
        %Time Matrix Initiation
        xlsPath = fullfile(xlsDir,xlsFolders{1},[xlsFolders{1} '_data.xlsx']);
        Type(p).Data(d).time = xlsread(xlsPath, 'A2:A314'); %Time matrix - not dependent on fly
        Type(p).Data(d).time = Type(p).Data(d).time(1:length(damXYZ)); %Creates time matrix that is the same length as the XYZ points
        
        
        %Time of Contact
        xlsPath2 = fullfile(xlsDir2,xlsFolders2{1},[xlsFolders2{1} '_data.xlsx']); %sets frame of contact/escape for each occurence
        tContactKir21Capture = xlsread(xlsPath2, 'A2:A14');
        tContactKir21Escape = xlsread(xlsPath2, 'A15:A20');
        tContactDLCapture = xlsread(xlsPath2, 'A21:A25');
        tContactDLEscape = xlsread(xlsPath2, 'A26:A34');
        
        
        %Pixels to mm conversion
        damXYZmm = damXYZ*25.4;
        flyXYZmm = flyXYZ*25.4;
        
        
        %Shortens time matrix to the point of contact/escape
        if p == 1,  Type(p).Data(d).t =  Type(p).Data(d).time(1:tContactKir21Capture(d));
        elseif p == 2,  Type(p).Data(d).t =  Type(p).Data(d).time(1:tContactKir21Escape(d));
        elseif p == 3, Type(p).Data(d).t = Type(p).Data(d).time(1:tContactDLCapture(d));
        elseif p == 4, Type(p).Data(d).t = Type(p).Data(d).time(1:tContactDLEscape(d));
        end
        
        
        %% STEP 2A: Butterworth Filter on each component axis
        
         acquisition_rate = 1000; % frames per sec recorded
        fNorm = 10 / (acquisition_rate/2);
        [b, a] = butter(2, fNorm);
        
        %Filtered XYZ of damselfly
        Type(p).Data(d).X_of_dam = filtfilt(b, a, damXYZmm(:,1));
        Type(p).Data(d).Y_of_dam = filtfilt(b, a, damXYZmm(:,2));
        Type(p).Data(d).Z_of_dam = filtfilt(b, a, damXYZmm(:,3));

        %Filtered XYZ of fly
        Type(p).Data(d).X_of_fly = filtfilt(b, a, flyXYZmm(:,1));
        Type(p).Data(d).Y_of_fly = filtfilt(b, a, flyXYZmm(:,2));
        Type(p).Data(d).Z_of_fly = filtfilt(b, a, flyXYZmm(:,3));
% 
% windowwidth = 15;
% b = ones(windowwidth,1)/windowwidth;
% Type(p).Data(d).X_of_dam = filter(b,1,damXYZmm(:,1));
% Type(p).Data(d).Y_of_dam = filter(b,1,damXYZmm(:,2));
% Type(p).Data(d).Z_of_dam = filter(b,1,damXYZmm(:,3));
% 
% Type(p).Data(d).X_of_dam(1:windowwidth) = NaN;
% Type(p).Data(d).Y_of_dam(1:windowwidth) = NaN;
% Type(p).Data(d).Z_of_dam(1:windowwidth) = NaN;

% Type(p).Data(d).X_of_fly = filter(kernal,1,flyXYZmm(:,1));
% Type(p).Data(d).Y_of_fly = filter(kernal,1,flyXYZmm(:,2));
% Type(p).Data(d).Z_of_fly = filter(kernal,1,flyXYZmm(:,3));
        
        %% STEP 2B: Example Plot of X, Y, and Z after Butterworth Filter
        
        %         ax1 = subplot(3,1,1); % top subplot
        %         ax2 = subplot(3,1,2); % middle subplot
        %         ax3 = subplot(3,1,3); % bottom subplot
        %
        %         plot(ax1, Type(p).Data(d).time, Type(p).Data(d).X_of_dam, Type(p).Data(d).time, damXYZmm(:,1))
        %
        %         plot(ax2, Type(p).Data(d).time, Type(p).Data(d).Y_of_dam, Type(p).Data(d).time, damXYZmm(:,2))
        %
        %         plot(ax3, Type(p).Data(d).time, Type(p).Data(d).Z_of_dam, Type(p).Data(d).time, damXYZmm(:,3))
        
        %% STEP 3: Distance Between Damselfly and Fly
        
        Type(p).Data(d).X_dist = Type(p).Data(d).X_of_dam/1000 - Type(p).Data(d).X_of_fly/1000;
        Type(p).Data(d).Y_dist = Type(p).Data(d).Y_of_dam/1000 - Type(p).Data(d).Y_of_fly/1000;
        Type(p).Data(d).Z_dist = Type(p).Data(d).Z_of_dam/1000 - Type(p).Data(d).Z_of_fly/1000;
        
        Type(p).Data(d).distXYZ = sqrt(Type(p).Data(d).X_dist.^2 + Type(p).Data(d).Y_dist.^2 + Type(p).Data(d).Z_dist.^2); %distance in m
        
        
        %% STEP 4: Translated Damselfly and Fly position relative to Fly Starting Position in m
        
        Type(p).Data(d).dam_trans_X = Type(p).Data(d).X_of_dam/1000 - ((Type(p).Data(d).tail(1) * 25.4)/1000);
        Type(p).Data(d).dam_trans_Y = Type(p).Data(d).Y_of_dam/1000 - ((Type(p).Data(d).tail(2) * 25.4)/1000);
        Type(p).Data(d).dam_trans_Z = Type(p).Data(d).Z_of_dam/1000 - ((Type(p).Data(d).tail(3) * 25.4)/1000);
        
        Type(p).Data(d).fly_trans_X = Type(p).Data(d).X_of_fly/1000 - ((Type(p).Data(d).tail(1) * 25.4)/1000);
        Type(p).Data(d).fly_trans_Y = Type(p).Data(d).Y_of_fly/1000 - ((Type(p).Data(d).tail(2) * 25.4)/1000);
        Type(p).Data(d).fly_trans_Z = Type(p).Data(d).Z_of_fly/1000 - ((Type(p).Data(d).tail(3) * 25.4)/1000);
        
        
        %% STEP 5: Azimuth and Elevation Calculations of Damselfly Rotated Such That Fly Orientation is Along the X-Axis
        
        % Translate fly axis points to tail at (0,0,0)
        Type(p).Data(d).head_trans = ((Type(p).Data(d).head * 25.4)/1000) - ((Type(p).Data(d).tail * 25.4)/1000); %in m
        Type(p).Data(d).tail_trans = ((Type(p).Data(d).tail * 25.4)/1000) - ((Type(p).Data(d).tail * 25.4)/1000); %in m
        
        %Solve for y-axis rotation angle
        z = Type(p).Data(d).head_trans(3);
        x1 = Type(p).Data(d).head_trans(1);
        theta = atan2(z,x1); %rotation by this angle places the fly orientation vector in the xy plane (z = 0)
        
        %Apply y-axis rotation angle to the damselfly points and fly axis
        damXtemp = Type(p).Data(d).dam_trans_X.*cos(theta) + Type(p).Data(d).dam_trans_Z.*sin(theta);
        damYtemp = Type(p).Data(d).dam_trans_Y;
        damZtemp = Type(p).Data(d).dam_trans_Z.*cos(theta) - Type(p).Data(d).dam_trans_X.*sin(theta);
        
        fly1Xtemp = Type(p).Data(d).fly_trans_X.*cos(theta) + Type(p).Data(d).fly_trans_Z.*sin(theta);
        fly1Ytemp = Type(p).Data(d).fly_trans_Y;
        fly1Ztemp = Type(p).Data(d).fly_trans_Z.*cos(theta) - Type(p).Data(d).fly_trans_X.*sin(theta);
        
        flyXtemp = Type(p).Data(d).head_trans(1).*cos(theta) + Type(p).Data(d).head_trans(3).*sin(theta);
        flyYtemp = Type(p).Data(d).head_trans(2);
        flyZtemp = Type(p).Data(d).head_trans(3).*cos(theta) - Type(p).Data(d).head_trans(1).*sin(theta);
        
        %Solve for z-axis rotation angle
        y = flyYtemp;
        x2 = flyXtemp;
        rho = -atan2(y,x2); %rotation by this angle places the fly orientation vector along the x-axis (y = 0, z = 0)
        
        %Apply z-axis rotation angle to the damselfly points
        Type(p).Data(d).dam_rot_X = damXtemp.*cos(rho) - damYtemp.*sin(rho);
        Type(p).Data(d).dam_rot_Y = damYtemp.*cos(rho) + damXtemp.*sin(rho);
        Type(p).Data(d).dam_rot_Z = damZtemp;
        
        Type(p).Data(d).fly_rot_X = fly1Xtemp.*cos(rho) - fly1Ytemp.*sin(rho);
        Type(p).Data(d).fly_rot_Y = fly1Ytemp.*cos(rho) + fly1Xtemp.*sin(rho);
        Type(p).Data(d).fly_rot_Z = fly1Ztemp;
        
        Type(p).Data(d).flyX = flyXtemp.*cos(rho) - flyYtemp.*sin(rho);
        Type(p).Data(d).flyY = flyYtemp.*cos(rho) + flyXtemp.*sin(rho);
        Type(p).Data(d).flyZ = flyZtemp;
        
        
        %Translate coordinates so fly location is (0,0,0)
        Type(p).Data(d).dam_rot_X = Type(p).Data(d).dam_rot_X - Type(p).Data(d).fly_rot_X;
        Type(p).Data(d).dam_rot_Y = Type(p).Data(d).dam_rot_Y - Type(p).Data(d).fly_rot_Y;
        Type(p).Data(d).dam_rot_Z = Type(p).Data(d).dam_rot_Z - Type(p).Data(d).fly_rot_Z;
        
        %The following requires all trajectories to start between 0 and 90 degree azimuth
        %After the start, the trajectories can be outside of this range
        %         if Type(p).Data(d).dam_rot_X(1)< 0
        %             Type(p).Data(d).dam_rot_X = -1.*Type(p).Data(d).dam_rot_X;
        %         end
        %         if Type(p).Data(d).dam_rot_Z(1)< 0
        %             Type(p).Data(d).dam_rot_Z = -1.*Type(p).Data(d).dam_rot_Z;
        %         end
        %
        %The following requires all data to be between 0 and 90 degree azimuth
        %This is for the duration of the trajectory
        Type(p).Data(d).dam_rot_X2 = Type(p).Data(d).dam_rot_X;
        Type(p).Data(d).dam_rot_Y2 = Type(p).Data(d).dam_rot_Y;
        Type(p).Data(d).dam_rot_Z2 = Type(p).Data(d).dam_rot_Z;
        
        Type(p).Data(d).dam_rot_X2(Type(p).Data(d).dam_rot_X2<0) = -1.*Type(p).Data(d).dam_rot_X2(Type(p).Data(d).dam_rot_X2<0);
        %         Type(p).Data(d).dam_rot_Z2(Type(p).Data(d).dam_rot_Z2<0) = -1.*Type(p).Data(d).dam_rot_Z2(Type(p).Data(d).dam_rot_Z2<0);
        
        
        [Type(p).Data(d).damAzi, Type(p).Data(d).damEle, ~] = cart2sph(Type(p).Data(d).dam_rot_X2, Type(p).Data(d).dam_rot_Y2, Type(p).Data(d).dam_rot_Z2); % Y and Z are swapped because fly is on the wall in the x-z plane
        
        Type(p).Data(d).damAzi = Type(p).Data(d).damAzi(1:length(Type(p).Data(d).t));
        Type(p).Data(d).damEle = Type(p).Data(d).damEle(1:length(Type(p).Data(d).t));
        
        avgAzimuth{p}(d) = circ_mean(Type(p).Data(d).damAzi(1:length(Type(p).Data(d).t))); % mean azimuth angle of damselfly up until time of contact
        
        avgElevation{p}(d) = circ_mean(Type(p).Data(d).damEle(1:length(Type(p).Data(d).t))); % mean elevation angle of damselfly up until time of contact
        
        avgAzimuth{p}(d) = rad2deg(avgAzimuth{p}(d)); % in deg
        
        avgElevation{p}(d) = rad2deg(avgElevation{p}(d)); % in deg
        
        %% STEP 5B: Plot rotated and non-rotated trajectories
        
        %The following is a check to make sure that the calculated distance between damselfly and fly initial position in unchanged between the rotated and unrotated points
        %         Type(p).Data(d).distXYZ_trans = sqrt(Type(p).Data(d).dam_trans_X.^2 + Type(p).Data(d).dam_trans_Y.^2 + Type(p).Data(d).dam_trans_Z.^2); %distance in m
        %         Type(p).Data(d).distXYZ_rot = sqrt(Type(p).Data(d).dam_rot_X.^2 + Type(p).Data(d).dam_rot_Y.^2 + Type(p).Data(d).dam_rot_Z.^2); %distance in m
        
        %Check to see if rotated and unrotated trajectories have similar shapes
        
        if d == numel(csvFolders) %plot all trajectories once the last is calculated in a group
            figure
            for i = 1:d
                plot3(Type(p).Data(i).dam_trans_X(1:length(Type(p).Data(i).t)),Type(p).Data(i).dam_trans_Y(1:length(Type(p).Data(i).t)),Type(p).Data(i).dam_trans_Z(1:length(Type(p).Data(i).t)))
                hold on
            end
            plot3(0,0,0,'o')
            
            
            if p == 1, title('Kir2.1(GF Silenced) Damselfly Raw Trajectories - Capture');
            elseif p == 2, title('Kir2.1 (GF Silenced) Damselfly Raw Trajectories - Escape');
            elseif p == 3, title('DL (GF Wild Type) Damselfly Raw Trajectories - Capture');
            elseif p == 4, title('DL (GF Wild Type) Damselfly Raw Trajectories - Escape');
            end
            
            xlabel('X (m)')
            ylabel('Y (m)')
            zlabel('Z (m)')
            
            figure
            for i = 1:d
                plot(Type(p).Data(i).dam_trans_Y(1:length(Type(p).Data(i).t)),Type(p).Data(i).dam_trans_Z(1:length(Type(p).Data(i).t)))
                hold on
            end
            plot3(0,0,0,'o')
            
            
            if p == 1, title('Kir2.1(GF Silenced) Damselfly Raw Trajectories - Capture');
            elseif p == 2, title('Kir2.1 (GF Silenced) Damselfly Raw Trajectories - Escape');
            elseif p == 3, title('DL (GF Wild Type) Damselfly Raw Trajectories - Capture');
            elseif p == 4, title('DL (GF Wild Type) Damselfly Raw Trajectories - Escape');
            end
            set(gca, 'XDir','reverse')
            
            xlabel('Y (m)')
            ylabel('Z (m)')
            
            figure
            for i = 1:d
                plot3(Type(p).Data(i).dam_rot_X(1:length(Type(p).Data(i).t)),Type(p).Data(i).dam_rot_Y(1:length(Type(p).Data(i).t)),Type(p).Data(i).dam_rot_Z(1:length(Type(p).Data(i).t)))
                hold on
            end
            plot3(0,0,0,'o')
            if p == 1, title('Kir2.1(GF Silenced) Damselfly Rotated Trajectories - Capture');
            elseif p == 2, title('Kir2.1 (GF Silenced) Damselfly Rotated Trajectories - Escape');
            elseif p == 3, title('DL (GF Wild Type) Damselfly Rotated Trajectories - Capture');
            elseif p == 4, title('DL (GF Wild Type) Damselfly Rotated Trajectories - Escape');
            end
            xlabel('X (m)')
            ylabel('Y (m)')
            zlabel('Z (m)')
        end
        
        
        
        %% STEP 6: Sgolay Differentiation of X, Y, and Z (Calculation of Smoothed Velocity and Acceleration)
        
        window_leg = 5;
        window = window_leg*2+1;
        
        headSize = 0.0014809; %radius of damselfly head in m
%         
%                 [b,g] = sgolay(3,window);
%         
%                 dx = zeros(length((Type(p).Data(d).X_of_dam/1000-Type(p).Data(d).X_of_fly(1)/1000)),3);
%                 dy = zeros(length((Type(p).Data(d).Y_of_dam/1000-Type(p).Data(d).Y_of_fly(1)/1000)),3);
%                 dz = zeros(length((Type(p).Data(d).Z_of_dam/1000-Type(p).Data(d).Z_of_fly(1)/1000)),3);
%                 dt = 1/acquisition_rate;
%                 for ind = 0:2
%                     dx(:,ind+1) = conv((Type(p).Data(d).X_of_dam/1000-Type(p).Data(d).X_of_fly(1)/1000), factorial(ind)/(-dt)^ind * g(:,ind+1), 'same');
%                     dy(:,ind+1) = conv((Type(p).Data(d).Y_of_dam/1000-Type(p).Data(d).Y_of_fly(1)/1000), factorial(ind)/(-dt)^ind * g(:,ind+1), 'same');
%                     dz(:,ind+1) = conv((Type(p).Data(d).Z_of_dam/1000-Type(p).Data(d).Z_of_fly(1)/1000), factorial(ind)/(-dt)^ind * g(:,ind+1), 'same');
%                 end
%         
%                 distVec_X = dx(:,1);
%                 speedVec_X = dx(:,2);
%                 accelVec_X = dx(:,3);
%         
%                 distVec_Y = dy(:,1);
%                 speedVec_Y = dy(:,2);
%                 accelVec_Y = dy(:,3);
%         
%                 distVec_Z = dz(:,1);
%                 speedVec_Z = dz(:,2);
%                 accelVec_Z = dz(:,3);
        
        
%         distVec_X = (Type(p).Data(d).X_of_dam/1000-Type(p).Data(d).X_of_fly(1)/1000);
%         distVec_Y = (Type(p).Data(d).Y_of_dam/1000-Type(p).Data(d).Y_of_fly(1)/1000);
%         distVec_Z = (Type(p).Data(d).Z_of_dam/1000-Type(p).Data(d).Z_of_fly(1)/1000);
%         
%         speedVec_X = zeros(length(Type(p).Data(d).X_of_dam),1);
%         speedVec_Y = zeros(length(Type(p).Data(d).X_of_dam),1);
%         speedVec_Z = zeros(length(Type(p).Data(d).X_of_dam),1);
%         
%         accelVec_X = zeros(length(Type(p).Data(d).X_of_dam),1);
%         accelVec_Y = zeros(length(Type(p).Data(d).X_of_dam),1);
%         accelVec_Z = zeros(length(Type(p).Data(d).X_of_dam),1);
%         
%         for i = (2:length(Type(p).Data(d).X_of_dam))
%             speedVec_X(i) = abs(distVec_X(i)-distVec_X(i-1))*acquisition_rate;
%             speedVec_Y(i) = abs(distVec_Y(i)-distVec_Y(i-1))*acquisition_rate;
%             speedVec_Z(i) = abs(distVec_Z(i)-distVec_Z(i-1))*acquisition_rate;
%         end
%         
%         for i = 2:length(Type(p).Data(d).X_of_dam)
%             accelVec_X(i) = speedVec_X(i)-speedVec_X(i-1)*acquisition_rate;
%             accelVec_Y(i) = speedVec_Y(i)-speedVec_Y(i-1)*acquisition_rate;
%             accelVec_Z(i) = speedVec_Z(i)-speedVec_Z(i-1)*acquisition_rate;
%         end
%         
        
%                 [distVec_X,speedVec_X,accelVec_X] = golayDifferentiate((Type(p).Data(d).X_of_dam/1000-Type(p).Data(d).X_of_fly(1)/1000),window);
%                 speedVec_X = speedVec_X*acquisition_rate;% m/s
%                 accelVec_X = accelVec_X*(acquisition_rate.^2);% m/s^2
%                 distVec_X(1:window_leg) = NaN;
%                 speedVec_X(1:window_leg) = NaN;
%                 accelVec_X(1:window_leg) = NaN;
%         
%                 [distVec_Y,speedVec_Y,accelVec_Y] = golayDifferentiate((Type(p).Data(d).Y_of_dam/1000-Type(p).Data(d).Y_of_fly(1)/1000),window);
%                 speedVec_Y = speedVec_Y*acquisition_rate;% m/s
%                 accelVec_Y = accelVec_Y*(acquisition_rate.^2);% m/s^2
%                 distVec_Y(1:window_leg) = NaN;
%                 speedVec_Y(1:window_leg) = NaN;
%                 accelVec_Y(1:window_leg) = NaN;
%         
%                 [distVec_Z,speedVec_Z,accelVec_Z] = golayDifferentiate((Type(p).Data(d).Z_of_dam/1000-Type(p).Data(d).Z_of_fly(1)/1000),window);
%                 speedVec_Z = speedVec_Z*acquisition_rate;% m/s
%                 accelVec_Z = accelVec_Z*(acquisition_rate.^2);% m/s^2
%                 distVec_Z(1:window_leg) = NaN;
%                 speedVec_Z(1:window_leg) = NaN;
%                 accelVec_Z(1:window_leg) = NaN;

distVec = sqrt((Type(p).Data(d).X_of_dam/1000-Type(p).Data(d).X_of_fly(1)/1000).^2 + (Type(p).Data(d).Y_of_dam/1000-Type(p).Data(d).Y_of_fly(1)/1000).^2 + (Type(p).Data(d).Z_of_dam/1000-Type(p).Data(d).Z_of_fly(1)/1000).^2);

distTrav = zeros(length(distVec),1);
for i = 2:length(distVec)
distTrav(i) = distTrav(i-1) + abs(distVec(i)-distVec(i-1));
end
                [distVec_Dams,~,~] = golayDifferentiate(distVec,window);
                [~,speedVec_Dams,accelVec_Dams] = golayDifferentiate(distTrav,window);
                speedVec_Dams = speedVec_Dams*acquisition_rate;% m/s
                accelVec_Dams = accelVec_Dams*(acquisition_rate.^2);% m/s^2
                distVec_Dams(1:window_leg) = NaN;
                speedVec_Dams(1:window_leg) = NaN;
                accelVec_Dams(1:window_leg) = NaN;
                
        
        
        %% STEP 7a: Combining X, Y, and Z Coordinates for distance
        
        %Type(p).Data(d).distVec_Dams = sqrt(distVec_X.^2 + distVec_Y.^2 + distVec_Z.^2);
        Type(p).Data(d).distVec_Dams = distVec_Dams(1:length(Type(p).Data(d).t));
% Type(p).Data(d).distVec_Dams = distTrav(1:length(Type(p).Data(d).t));
        
        %% STEP 7: Combining X, Y, and Z Coordiantes for Velocity - Peak Speed, Size to Speed Ratios, Average Speed
        
        %Type(p).Data(d).speedVec_Dams = sqrt(speedVec_X.^2 + speedVec_Y.^2 + speedVec_Z.^2);
        Type(p).Data(d).speedVec_Dams = speedVec_Dams(1:length(Type(p).Data(d).t));
        
        maxSpeed{p}(d) = nanmax(Type(p).Data(d).speedVec_Dams(1:length(Type(p).Data(d).t))); % max speed until TOC in m/s
        sizeSpeed{p}(d) = (headSize/maxSpeed{p}(d)) * 1000; % in ms
        initialSizeSpeed{p}(d) = (headSize/ Type(p).Data(d).speedVec_Dams(8)) * 1000; % in ms
        
        avgSpeed{p}(d) = nanmean(Type(p).Data(d).speedVec_Dams(1:length(Type(p).Data(d).t))); % average speed until TOC
        
        avgSizeSpeed{p}(d) = (headSize/avgSpeed{p}(d)) * 1000; % in ms
        
        Type(p).Data(d).sizeSpeedtable = headSize./Type(p).Data(d).speedVec_Dams(1:length(Type(p).Data(d).t)) * 1000; %Size to speed ratio for each time step
        
        
        %% STEP 8: Combining X, Y, and Z coordinates for acceleration
        
        %Type(p).Data(d).accelVec_Dams = sqrt(accelVec_X.^2 + accelVec_Y.^2 + accelVec_Z.^2);
        Type(p).Data(d).accelVec_Dams = accelVec_Dams(1:length(Type(p).Data(d).t));
        
        peakAcc{p}(d) = nanmax(Type(p).Data(d).accelVec_Dams); % max acceleration until TOC in m/s^2
        
        accAvg{p}(d) = nanmean(Type(p).Data(d).accelVec_Dams(1:length(Type(p).Data(d).t))); % average acceleration until TOC
        
        % If there are more than 68 time steps (first 8 are NaN), accAvg2 is the average acceleration in the last 60 time steps
        if length(Type(p).Data(d).t) <= 68, accAvg2{p}(d) = nanmean(Type(p).Data(d).accelVec_Dams(1:length(Type(p).Data(d).t)));
        elseif length(Type(p).Data(d).t) > 68, accAvg2{p}(d) = nanmean(Type(p).Data(d).accelVec_Dams(length(Type(p).Data(d).t) - 60:length(Type(p).Data(d).t)));
        end
        
        %Same as above, but for median acceleration
        if length(Type(p).Data(d).t) <= 68, accMed2{p}(d) = nanmedian(Type(p).Data(d).accelVec_Dams(1:length(Type(p).Data(d).t)));
        elseif length(Type(p).Data(d).t) > 68, accMed2{p}(d) = nanmedian(Type(p).Data(d).accelVec_Dams(length(Type(p).Data(d).t) - 60:length(Type(p).Data(d).t)));
        end
        
        
        %% STEP 9: Tier Plots (distance, speed, acceleration)
        
        %         Type(p).Data(d).distVec_Dams = sqrt((distVec_X/1000).^2 + (distVec_Y/1000).^2 + (distVec_Z/1000).^2);
        %
        %         figure(10)
        %         ax1 = subplot(3,1,1); % top subplot
        %         ax2 = subplot(3,1,2); % middle subplot
        %         ax3 = subplot(3,1,3); % bottom subplot
        %
        %         plot(ax1,t,Type(p).Data(d).distVec_Dams(1:length(t)))
        %         title(ax1,'Position');
        %         xlabel(ax1, 'Time (s)');
        %         ylabel(ax1, 'Position (m)');
        %
        %         plot(ax2,t,Type(p).Data(d).speedVec_Dams(1:length(t)))
        %         title(ax2,'Speed');
        %         xlabel(ax2, 'Time (s)');
        %         ylabel(ax2, 'Speed (m/s)');
        %         refline(ax2, [0 maxSpeed{p}(d)]);
        %         text(ax2, 0.02, maxSpeed{p}(d) - 0.02, 'Peak Speed (m/s)', 'fontsize', 10);
        %
        %         plot(ax3,t,Type(p).Data(d).accelVec_Dams(1:length(t)))
        %         title(ax3,'Acceleration');
        %         xlabel(ax3, 'Time (s)');
        %         ylabel(ax3, 'Acceleration (m/s^2)');
        
        %         export_fig('Tier Plots.pdf', '-append')
        
        
        %% STEP 10: Tier Plots for Azimuth and Elevation
        %
        if d == numel(csvFolders) %plot all trajectories once the last is calculated in a group
            figure
            for i = 1:d
                plot(Type(p).Data(i).t,rad2deg(Type(p).Data(i).damAzi(1:length(Type(p).Data(i).t))))
                if p == 1, title('Kir2.1(GF Silenced) Damselfly Azimuth - Capture');
                elseif p == 2, title('Kir2.1 (GF Silenced) Damselfly Azimuth - Escape');
                elseif p == 3, title('DL (GF Wild Type) Damselfly Azimuth - Capture');
                elseif p == 4, title('DL (GF Wild Type) Damselfly Azimuth - Escape');
                end
                xlabel('Time (s)')
                ylabel('Azimuth (ï¿½)')
                hold on
            end
            
            figure
            for i = 1:d
                plot(Type(p).Data(i).t, rad2deg(Type(p).Data(i).damEle(1:length(Type(p).Data(i).t))))
                if p == 1, title('Kir2.1(GF Silenced) Damselfly Elevation - Capture');
                elseif p == 2, title('Kir2.1 (GF Silenced) Damselfly Elevation - Escape');
                elseif p == 3, title('DL (GF Wild Type) Damselfly Elevation - Capture');
                elseif p == 4, title('DL (GF Wild Type) Damselfly Elevation - Escape');
                end
                xlabel('Time (s)')
                ylabel('Elevation (ï¿½)')
                hold on
            end
        end
        %
        %         export_fig('Azimuth_Elevation_Tier Plots.pdf', '-append')
        
        
        %% STEP 11: OTHER KINEMATICS
        
        
        
        tvVisualAng{p}(d) = 2 * atand((headSize)/(Type(p).Data(d).distXYZ(2) - (Type(p).Data(d).speedVec_Dams(8))*((Type(p).Data(d).t(length(Type(p).Data(d).t))) - (0.5)*(accAvg{p}(d))*((Type(p).Data(d).t(length(Type(p).Data(d).t)))^2))));
        
        
        %Angular Size of Damselfly
        
        Type(p).Data(d).angSizetable = 2 * atand(headSize./Type(p).Data(d).distVec_Dams(8:length(Type(p).Data(d).t))); %Angular size of Damselfly over time
        
        %angSize{p}(d) = Type(p).Data(d).angSizetable(end); %Angular size of Damselfly at TOC
        
        
        timeOfCont{p}(d) = Type(p).Data(d).t(end); %Time of Contact
        
        distTOC{p}(d) = Type(p).Data(d).distVec_Dams(length(Type(p).Data(d).t)); % distance to fly at TOC
        
        
        angSize{p}(d) =  2 * atand(headSize./Type(p).Data(d).distVec_Dams(length(Type(p).Data(d).t)));
        eventdur{p}(d) = length(Type(p).Data(d).t); %in ms
        
    end
    
    
end

%% STEP 12: Plot all rotated trajectories together
%linecolor = [[1 .2 .2] ;[.2 .2 1] ;[.6 0 0] ;[0 0 .6]];
linecolor = [[1 0 0] ;[0 0 1] ;[1 .5 0] ;[0 1 0]];
figure('Position',[100 100 600 600])
for p = 1:4
    for i = 1:length(Type(p).Data)
        %         if Type(p).Data(i).dam_rot_X(1)< 0
        %             Type(p).Data(i).dam_rot_X = -1.*Type(p).Data(i).dam_rot_X;
        %         end
        %         if Type(p).Data(i).dam_rot_Z(1)< 0
        %             Type(p).Data(i).dam_rot_Z = -1.*Type(p).Data(i).dam_rot_Z;
        %         end
        if p == 1; a = plot3(Type(p).Data(i).dam_rot_X(1:length(Type(p).Data(i).t))*1000,Type(p).Data(i).dam_rot_Y(1:length(Type(p).Data(i).t))*1000,Type(p).Data(i).dam_rot_Z(1:length(Type(p).Data(i).t))*1000,'Color',linecolor(p,:));
        elseif p ==2; b = plot3(Type(p).Data(i).dam_rot_X(1:length(Type(p).Data(i).t))*1000,Type(p).Data(i).dam_rot_Y(1:length(Type(p).Data(i).t))*1000,Type(p).Data(i).dam_rot_Z(1:length(Type(p).Data(i).t))*1000,'Color',linecolor(p,:));
        elseif p ==3; c = plot3(Type(p).Data(i).dam_rot_X(1:length(Type(p).Data(i).t))*1000,Type(p).Data(i).dam_rot_Y(1:length(Type(p).Data(i).t))*1000,Type(p).Data(i).dam_rot_Z(1:length(Type(p).Data(i).t))*1000,'Color',linecolor(p,:));
        else; d = plot3(Type(p).Data(i).dam_rot_X(1:length(Type(p).Data(i).t))*1000,Type(p).Data(i).dam_rot_Y(1:length(Type(p).Data(i).t))*1000,Type(p).Data(i).dam_rot_Z(1:length(Type(p).Data(i).t))*1000,'Color',linecolor(p,:));
        end
        hold on
        plot3(Type(p).Data(i).dam_rot_X(1)*1000,Type(p).Data(i).dam_rot_Y(1)*1000,Type(p).Data(i).dam_rot_Z(1)*1000,'o','MarkerEdgeColor',linecolor(p,:),'MarkerFaceColor',linecolor(p,:))
    end
end
e = plot3(0,0,0,'kd','MarkerFaceColor','k');
% title('Damselfly Trajectories');
[~, hobj, ~, ~] = legend([c(1) a(1) d(1) b(1) e(1)], 'GF1>+ - Capture','GF1>Kir2.1 - Capture','GF1>+ - Escape', 'GF1>Kir2.1 - Escape','Fly Location');
set(hobj,'LineWidth',1.5);
xlabel('X (mm)','FontWeight', 'bold', 'FontSize', 14)
ylabel('Y (mm)','FontWeight', 'bold', 'FontSize', 14)
zlabel('Z (mm)','FontWeight', 'bold', 'FontSize', 14)
ax = gca;
ax.FontSize = 14;
pbaspect([1 1 1])
set(gca, 'box', 'off')

h = rotate3d;
set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
set(gcf, 'ResizeFcn', @align_axislabel)
align_axislabel([], gca)
axislabel_translation_slider;

%% STEP 12B: Plot seperate axes for rotated trajectories
figure('Position',[100 100 600 500])
%figure
for p = 1:4
    for i = 1:length(Type(p).Data)
        plot(Type(p).Data(i).dam_rot_Y(1:length(Type(p).Data(i).t))*1000,Type(p).Data(i).dam_rot_X(1:length(Type(p).Data(i).t))*1000,'Color',linecolor(p,:));
        dataoutx = Type(p).Data(i).dam_rot_Y(1:length(Type(p).Data(i).t))*1000;
        dataoutz = Type(p).Data(i).dam_rot_X(1:length(Type(p).Data(i).t))*1000;
        dataout = [dataoutx dataoutz];
        hold on
        plot(Type(p).Data(i).dam_rot_Y(1)*1000,Type(p).Data(i).dam_rot_X(1)*1000,'o','MarkerEdgeColor',linecolor(p,:),'MarkerFaceColor',linecolor(p,:))
    end
    a=0;
end
e = plot(0,0,'kd','MarkerFaceColor','k','MarkerSize',9);
%title('Damselfly Trajectories');
xlabel('X (mm)','FontWeight', 'bold', 'FontSize', 14)
ylabel('Y (mm)','FontWeight', 'bold', 'FontSize', 14)
zlabel('Y (mm)','FontWeight', 'bold', 'FontSize', 14)
ylim([-20 30])
xlim([0 61])
ax = gca;
ax.FontSize = 14;
axis equal
%pbaspect([1 1 1])
set(gca, 'box', 'off')

%% STEP 13: Speed Plots Comparing Capture vs Escape

%Setup speed variables for plotting
speedp1 = padcat(flipud(Type(1).Data(1).speedVec_Dams), flipud(Type(1).Data(2).speedVec_Dams), flipud(Type(1).Data(3).speedVec_Dams), flipud(Type(1).Data(4).speedVec_Dams), flipud(Type(1).Data(5).speedVec_Dams)...
    ,flipud(Type(1).Data(6).speedVec_Dams), flipud(Type(1).Data(7).speedVec_Dams), flipud(Type(1).Data(8).speedVec_Dams), flipud(Type(1).Data(6).speedVec_Dams), flipud(Type(1).Data(10).speedVec_Dams), flipud(Type(1).Data(11).speedVec_Dams)...
    ,flipud(Type(1).Data(12).speedVec_Dams), flipud(Type(1).Data(13).speedVec_Dams));
speedp1 = flipud(speedp1); %concat speed matricies
numnansp1 = sum(~isnan(speedp1),2); %number of traces averaged as a function of time
meanp1Vec = nanmean(speedp1, 2); %mean speed as a function of time
stdp1Vec = nanstd(speedp1, 0, 2); %standard deviation of speed

speedp2 = padcat(flipud(Type(2).Data(1).speedVec_Dams), flipud(Type(2).Data(2).speedVec_Dams), flipud(Type(2).Data(3).speedVec_Dams), flipud(Type(2).Data(4).speedVec_Dams), flipud(Type(2).Data(5).speedVec_Dams)...
    , flipud(Type(2).Data(6).speedVec_Dams));
speedp2 = flipud(speedp2);
numnansp2 = sum(~isnan(speedp2),2);
meanp2Vec = nanmean(speedp2, 2);
stdp2Vec = nanstd(speedp2, 0, 2);

speedp3 = padcat(flipud(Type(3).Data(1).speedVec_Dams), flipud(Type(3).Data(2).speedVec_Dams), flipud(Type(3).Data(3).speedVec_Dams), flipud(Type(3).Data(4).speedVec_Dams), flipud(Type(3).Data(5).speedVec_Dams));
speedp3 = flipud(speedp3);
numnansp3 = sum(~isnan(speedp3),2);
meanp3Vec = nanmean(speedp3, 2);
stdp3Vec = nanstd(speedp3, 0, 2);

speedp4 = padcat(flipud(Type(4).Data(1).speedVec_Dams), flipud(Type(4).Data(2).speedVec_Dams), flipud(Type(4).Data(3).speedVec_Dams), flipud(Type(4).Data(4).speedVec_Dams), flipud(Type(4).Data(5).speedVec_Dams)...
    ,flipud(Type(4).Data(6).speedVec_Dams), flipud(Type(4).Data(7).speedVec_Dams) );
speedp4 = flipud(speedp4);
numnansp4 = sum(~isnan(speedp4),2);
meanp4Vec = nanmean(speedp4, 2);
stdp4Vec = nanstd(speedp4, 0, 2);


linecolor = [[1 0 0] ;[0 0 1] ;[1 .5 0] ;[0 1 0]];
% Plot DL Speeds - Capture and Escape
figure('Position',[500 300 700 635])
s1 = subplot(2, 1, 1);
for p = 1:length(Type)
    
    hold on
    for d = 1:length(Type(p).Data)
        
        %                     plot((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000, meanp3Vec(1:length(Type(3).Data(1).t)), 'r', 'LineWidth', 2); %mean
        %                     plot((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000, (meanp3Vec(1:length(Type(3).Data(1).t)) - stdp3Vec(1:length(Type(3).Data(1).t))), '--r'); %lower standard dev line
        %                     plot((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000, (meanp3Vec(1:length(Type(3).Data(1).t)) + stdp3Vec(1:length(Type(3).Data(1).t))), '--r') %upper standard dev line
        %
        %                     plot((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000, meanp4Vec(1:length(Type(4).Data(6).t)), 'b', 'LineWidth', 2);
        %                     plot((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000, (meanp4Vec(1:length(Type(4).Data(6).t)) - stdp4Vec(1:length(Type(4).Data(6).t))), '--b');
        %                     plot((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000, (meanp4Vec(1:length(Type(4).Data(6).t)) + stdp4Vec(1:length(Type(4).Data(6).t))), '--b')
        %
        % Uncomment for raw traces
                            if p == 3,  one = plot((Type(p).Data(d).t - Type(p).Data(d).t(end)) * 1000,Type(p).Data(d).speedVec_Dams(1:length(Type(p).Data(d).t)), 'Color',linecolor(p,:));
                            elseif p == 4, two = plot((Type(p).Data(d).t - Type(p).Data(d).t(end)) * 1000,Type(p).Data(d).speedVec_Dams(1:length(Type(p).Data(d).t)), 'Color',linecolor(p,:));
                            end
                            timedata = (Type(p).Data(d).t - Type(p).Data(d).t(end)) * 1000;
                            distdata = Type(p).Data(d).damEle(1:length(Type(p).Data(d).t));
                            combined = [timedata distdata];
                            a = 0;
    end
end

three = shadedErrorBar(((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000), meanp3Vec(1:length(Type(3).Data(1).t)), stdp3Vec(1:length(Type(3).Data(1).t)),'lineprops',{'Color',linecolor(3,:),'LineWidth',2},'transparent', 1); %shaded area for standard dev
hold on
four = shadedErrorBar(((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000), meanp4Vec(1:length(Type(4).Data(6).t)), stdp4Vec(1:length(Type(4).Data(6).t)), 'lineprops', {'Color',linecolor(4,:),'LineWidth',2}, 'transparent', 1);


title('GF1 > +','FontSize',14)
%xlabel('Time (ms)', 'FontWeight', 'bold','FontSize',12)
ylabel('Speed (m/s)', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[0,0.6],'Xlim',[-100 0])
set(s1,'Position',[.13 .5 .7750 .45])
hold off

s2 = subplot(2, 1, 2);
plot(((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000), numnansp3(1:length(Type(3).Data(1).t)), 'Color',linecolor(3,:)) %number of traces averaged
hold on, plot(((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000), numnansp4(1:length(Type(4).Data(6).t)), 'Color',linecolor(4,:))
set(gca,'Ylim', [0, 15], 'Xlim',[-100 0])

% legend([one two three four five six], 'Capture', 'Escape', 'Capture Mean', 'Escape Mean', 'Capture Std', 'Escape Std')
[~, hobj, ~, ~] = legend([three.mainLine,four.mainLine], 'Capture', 'Escape', 'Location', 'northwest', 'Orientation', 'horizontal');
set(hobj,'LineWidth',1.5);
xlabel('Time (ms)', 'FontWeight', 'bold','FontSize',14)
ylabel('n', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca, 'box', 'off')
set(s2,'Position',[.13 .33 .7750 .1])

% Plot Kir Speeds - Capture and Escape
figure('Position',[500 300 700 635])
s1 = subplot(2, 1, 1);
for p = 1:length(Type)
    
    hold on
    for d = 1:length(Type(p).Data)
        
        %
        %         if p == 1,  one = plot((Type(p).Data(d).t - Type(p).Data(d).t(end)) * 1000,Type(p).Data(d).speedVec_Dams(1:length(Type(p).Data(d).t)), 'Color',linecolor(p,:));
        %         elseif p == 2, two = plot((Type(p).Data(d).t - Type(p).Data(d).t(end)) * 1000,Type(p).Data(d).speedVec_Dams(1:length(Type(p).Data(d).t)), 'Color',linecolor(p,:));
        %         end
        
        %
        %                     three = plot((Type(1).Data(3).t - Type(1).Data(3).t(end))*1000, meanp1Vec(1:length(Type(1).Data(3).t)), 'r', 'LineWidth', 2);
        %
        %                     five = plot((Type(1).Data(3).t - Type(1).Data(3).t(end))*1000, (meanp1Vec(1:length(Type(1).Data(3).t)) - stdp1Vec(1:length(Type(1).Data(3).t))), '--r');
        %
        %                     plot((Type(1).Data(3).t - Type(1).Data(3).t(end))*1000, (meanp1Vec(1:length(Type(1).Data(3).t)) + stdp1Vec(1:length(Type(1).Data(3).t))), '--r')
        %
        %                     four = plot((Type(2).Data(5).t - Type(2).Data(5).t(end))*1000, meanp2Vec(1:length(Type(2).Data(5).t)), 'b', 'LineWidth', 2);
        %
        %                     six = plot((Type(2).Data(5).t - Type(2).Data(5).t(end))*1000, (meanp2Vec(1:length(Type(2).Data(5).t)) - stdp2Vec(1:length(Type(2).Data(5).t))), '--b');
        %
        %                     plot((Type(2).Data(5).t - Type(2).Data(5).t(end))*1000, (meanp2Vec(1:length(Type(2).Data(5).t)) + stdp2Vec(1:length(Type(2).Data(5).t))), '--b')
        
        
    end
    
end

three = shadedErrorBar((Type(1).Data(3).t - Type(1).Data(3).t(end))*1000, meanp1Vec(1:length(Type(1).Data(3).t)), stdp1Vec(1:length(Type(1).Data(3).t)), 'lineprops', {'Color',linecolor(1,:),'LineWidth',2}, 'transparent', 1);
hold on
four = shadedErrorBar((Type(2).Data(5).t - Type(2).Data(5).t(end))*1000, meanp2Vec(1:length(Type(2).Data(5).t)), stdp2Vec(1:length(Type(2).Data(5).t)), 'lineprops', {'Color',linecolor(2,:),'LineWidth',2}, 'transparent', 1);

title('GF1 > Kir2.1','FontSize',14)
%xlabel('Time (ms)', 'FontWeight', 'bold','FontSize',12)
ylabel('Speed (m/s)', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[0,0.5],'Xlim',[-100 0])
set(s1,'Position',[.13 .5 .7750 .45])
hold off

s2 = subplot(2, 1, 2);
plot(((Type(1).Data(3).t - Type(1).Data(3).t(end)) * 1000), numnansp1(1:length(Type(1).Data(3).t)), 'Color',linecolor(1,:))
hold on, plot(((Type(2).Data(5).t - Type(2).Data(5).t(end)) * 1000), numnansp2(1:length(Type(2).Data(5).t)), 'Color',linecolor(2,:))
set(gca,'Ylim', [0, 15], 'Xlim',[-100 0])
set(s2,'Position',[.13 .33 .7750 .1])

[~, hobj, ~, ~] = legend([three.mainLine,four.mainLine], 'Capture', 'Escape', 'Location', 'northwest', 'Orientation', 'horizontal');
set(hobj,'LineWidth',1.5);
xlabel('Time (ms)', 'FontWeight', 'bold','FontSize',14)
ylabel('n', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca, 'box', 'off')


%% STEP 14: Speed Plots Comparing DL vs Kir
figure
subplot(2, 1, 1);
for p = 1:length(Type)
    
    hold on
    for d = 1:length(Type(p).Data)
        
        if p == 1,  one = plot((Type(p).Data(d).t - Type(p).Data(d).t(end)) * 1000,Type(p).Data(d).speedVec_Dams(1:length(Type(p).Data(d).t)), 'g');
        elseif p == 3,  two = plot((Type(p).Data(d).t - Type(p).Data(d).t(end)) * 1000,Type(p).Data(d).speedVec_Dams(1:length(Type(p).Data(d).t)), 'b');
        end
        
        %         three = plot((Type(1).Data(3).t - Type(1).Data(3).t(end))*1000, meanp1Vec(1:length(Type(1).Data(3).t)), 'g', 'LineWidth', 2);
        %
        %         five = plot((Type(1).Data(3).t - Type(1).Data(3).t(end))*1000, (meanp1Vec(1:length(Type(1).Data(3).t)) - stdp1Vec(1:length(Type(1).Data(3).t))), '--g');
        %
        %         plot((Type(1).Data(3).t - Type(1).Data(3).t(end))*1000, (meanp1Vec(1:length(Type(1).Data(3).t)) + stdp1Vec(1:length(Type(1).Data(3).t))), '--g')
        %
        %         four = plot((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000, meanp3Vec(1:length(Type(3).Data(1).t)), 'b', 'LineWidth', 2);
        %
        %         six = plot((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000, (meanp3Vec(1:length(Type(3).Data(1).t)) - stdp3Vec(1:length(Type(3).Data(1).t))), '--b');
        %
        %         plot((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000, (meanp3Vec(1:length(Type(3).Data(1).t)) + stdp3Vec(1:length(Type(3).Data(1).t))), '--b')
        
        
        
    end
    
end

three = shadedErrorBar((Type(1).Data(3).t - Type(1).Data(3).t(end))*1000, meanp1Vec(1:length(Type(1).Data(3).t)), stdp1Vec(1:length(Type(1).Data(3).t)), 'lineprops', 'g', 'transparent', 1);
hold on
four = shadedErrorBar((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000, meanp3Vec(1:length(Type(3).Data(1).t)), stdp3Vec(1:length(Type(3).Data(1).t)), 'lineprops', 'b', 'transparent', 1);


hold on, title('Mean Capture Velocity Plots')
hold on, xlabel('Time(ms)', 'FontWeight', 'bold','FontSize',14)
hold on, ylabel('Speed (m/s)', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
%     hold on, set(gca,'Ylim',[0,0.5],'Xlim',[-100 0])
hold off

subplot(2, 1, 2);
plot(((Type(1).Data(3).t - Type(1).Data(3).t(end)) * 1000), numnansp1(1:length(Type(1).Data(3).t)), 'g')
hold on, plot(((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000), numnansp3(1:length(Type(3).Data(1).t)), 'b')
%     hold on, set(gca,'Ylim', [0, 14], 'Xlim',[-100 0])

% legend([one two three four five six], 'Kir2.1', 'DL', 'Kir2.1 Mean', 'DL Mean', 'Kir2.1 Std', 'DL Std')
legend([three.mainLine,four.mainLine], 'GF1 > Kir2.1', 'GF1 > +', 'Location', 'northwest', 'Orientation', 'horizontal')
figure

subplot(2, 1, 1);
for p = 1:length(Type)
    
    hold on
    for d = 1:length(Type(p).Data)
        
        if p == 2,  one = plot((Type(p).Data(d).t - Type(p).Data(d).t(end)) * 1000,Type(p).Data(d).speedVec_Dams(1:length(Type(p).Data(d).t)), 'g');
        elseif p == 4,  two = plot((Type(p).Data(d).t - Type(p).Data(d).t(end)) * 1000,Type(p).Data(d).speedVec_Dams(1:length(Type(p).Data(d).t)), 'b');
        end
        
        
        %         three = plot((Type(2).Data(5).t - Type(2).Data(5).t(end))*1000, meanp2Vec(1:length(Type(2).Data(5).t)), 'g', 'LineWidth', 2);
        %
        %         five = plot((Type(2).Data(5).t - Type(2).Data(5).t(end))*1000, (meanp2Vec(1:length(Type(2).Data(5).t)) - stdp2Vec(1:length(Type(2).Data(5).t))), '--g');
        %
        %         plot((Type(2).Data(5).t - Type(2).Data(5).t(end))*1000, (meanp2Vec(1:length(Type(2).Data(5).t)) + stdp2Vec(1:length(Type(2).Data(5).t))), '--g')
        %
        %
        %         four = plot((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000, meanp4Vec(1:length(Type(4).Data(6).t)), 'b', 'LineWidth', 2);
        %
        %         six = plot((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000, (meanp4Vec(1:length(Type(4).Data(6).t)) - stdp4Vec(1:length(Type(4).Data(6).t))), '--b');
        %
        %         plot((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000, (meanp4Vec(1:length(Type(4).Data(6).t)) + stdp4Vec(1:length(Type(4).Data(6).t))), '--b')
        
        
        
    end
    
end

three = shadedErrorBar((Type(2).Data(5).t - Type(2).Data(5).t(end))*1000, meanp2Vec(1:length(Type(2).Data(5).t)), stdp2Vec(1:length(Type(2).Data(5).t)), 'lineprops', 'g', 'transparent', 1);
hold on
four = shadedErrorBar((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000, meanp4Vec(1:length(Type(4).Data(6).t)), stdp4Vec(1:length(Type(4).Data(6).t)), 'lineprops', 'b', 'transparent', 1);

hold on, title('Mean Escape Velocity Plots')
hold on, xlabel('Time(ms)', 'FontWeight', 'bold','FontSize',14)
hold on, ylabel('Speed (m/s)', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
%     hold on, set(gca,'Ylim',[0,0.5],'Xlim',[-100 0])
hold off

subplot(2, 1, 2);
plot(((Type(2).Data(5).t - Type(2).Data(5).t(end)) * 1000), numnansp2(1:length(Type(2).Data(5).t)), 'g')
hold on, plot(((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000), numnansp4(1:length(Type(4).Data(6).t)), 'b')
%     hold on, set(gca,'Ylim', [0, 14], 'Xlim',[-100 0])

% legend([one two three four five six], 'Kir2.1', 'DL', 'Kir2.1 Mean', 'DL Mean', 'Kir2.1 Std', 'DL Std')
legend([three.mainLine,four.mainLine], 'GF1 > Kir2.1', 'GF1 > +', 'Location', 'northwest', 'Orientation', 'horizontal')
%% STEP 15: Speed Plots Combined
% Plot DL Speeds - Capture and Escape
figure('Position',[500 300 700 635])
s1 = subplot(2, 1, 1);
%
%                     plot((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000, meanp3Vec(1:length(Type(3).Data(1).t)),'Color',linecolor(3,:), 'LineWidth', 2); %mean
%                     plot((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000, (meanp3Vec(1:length(Type(3).Data(1).t)) - stdp3Vec(1:length(Type(3).Data(1).t))), '--','Color',linecolor(3,:)) %lower standard dev line
%                     plot((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000, (meanp3Vec(1:length(Type(3).Data(1).t)) + stdp3Vec(1:length(Type(3).Data(1).t))), '--','Color',linecolor(3,:)) %upper standard dev line
%
%                    plot((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000, meanp4Vec(1:length(Type(4).Data(6).t)), 'Color',linecolor(4,:), 'LineWidth', 2);
%                     plot((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000, (meanp4Vec(1:length(Type(4).Data(6).t)) - stdp4Vec(1:length(Type(4).Data(6).t))), '--','Color',linecolor(4,:))
%                     plot((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000, (meanp4Vec(1:length(Type(4).Data(6).t)) + stdp4Vec(1:length(Type(4).Data(6).t))), '--','Color',linecolor(4,:))
%
%                     plot((Type(1).Data(3).t - Type(1).Data(3).t(end))*1000, meanp1Vec(1:length(Type(1).Data(3).t)),'Color',linecolor(1,:), 'LineWidth', 2);
%                     plot((Type(1).Data(3).t - Type(1).Data(3).t(end))*1000, (meanp1Vec(1:length(Type(1).Data(3).t)) - stdp1Vec(1:length(Type(1).Data(3).t))), '--','Color',linecolor(1,:))
%                     plot((Type(1).Data(3).t - Type(1).Data(3).t(end))*1000, (meanp1Vec(1:length(Type(1).Data(3).t)) + stdp1Vec(1:length(Type(1).Data(3).t))), '--r','Color',linecolor(1,:))
%
%                    plot((Type(2).Data(5).t - Type(2).Data(5).t(end))*1000, meanp2Vec(1:length(Type(2).Data(5).t)),'Color',linecolor(2,:), 'LineWidth', 2);
%                     plot((Type(2).Data(5).t - Type(2).Data(5).t(end))*1000, (meanp2Vec(1:length(Type(2).Data(5).t)) - stdp2Vec(1:length(Type(2).Data(5).t))), '--','Color',linecolor(2,:))
%                     plot((Type(2).Data(5).t - Type(2).Data(5).t(end))*1000, (meanp2Vec(1:length(Type(2).Data(5).t)) + stdp2Vec(1:length(Type(2).Data(5).t))), '--b','Color',linecolor(2,:))
%
three = shadedErrorBar(((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000), meanp3Vec(1:length(Type(3).Data(1).t)), stdp3Vec(1:length(Type(3).Data(1).t)),'lineprops',{'Color',linecolor(3,:),'LineWidth',2},'transparent', 1); %shaded area for standard dev
hold on
four = shadedErrorBar(((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000), meanp4Vec(1:length(Type(4).Data(6).t)), stdp4Vec(1:length(Type(4).Data(6).t)), 'lineprops',{'Color', linecolor(4,:),'LineWidth',2}, 'transparent', 1);
hold on
one = shadedErrorBar((Type(1).Data(3).t - Type(1).Data(3).t(end))*1000, meanp1Vec(1:length(Type(1).Data(3).t)), stdp1Vec(1:length(Type(1).Data(3).t)), 'lineprops',{'Color', linecolor(1,:),'LineWidth',2}, 'transparent', 1);
hold on
two = shadedErrorBar((Type(2).Data(5).t - Type(2).Data(5).t(end))*1000, meanp2Vec(1:length(Type(2).Data(5).t)), stdp2Vec(1:length(Type(2).Data(5).t)), 'lineprops',{'Color', linecolor(2,:),'LineWidth',2}, 'transparent', 1);


%xlabel('Time (ms)', 'FontWeight', 'bold','FontSize',12)
ylabel('Speed (m/s)', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[0,0.5],'Xlim',[-100 0])
set(s1,'Position',[.13 .5 .7750 .45])
legend([one.mainLine,two.mainLine,three.mainLine,four.mainLine], 'GF1>Kir2.1 - Capture', 'GF1>Kir2.1 - Escape', 'GF1>+ - Capture', 'GF1>+ - Escape', 'Location', 'northwest', 'Orientation', 'horizontal');
set(legend, 'NumColumns' ,2)
set(gca, 'box', 'off')
hold off

s2 = subplot(2, 1, 2);
plot(((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000), numnansp3(1:length(Type(3).Data(1).t)), 'Color',linecolor(3,:)) %number of traces averaged
hold on, plot(((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000), numnansp4(1:length(Type(4).Data(6).t)), 'Color',linecolor(4,:))
hold on, plot(((Type(1).Data(3).t - Type(1).Data(3).t(end)) * 1000), numnansp1(1:length(Type(1).Data(3).t)), 'Color',linecolor(1,:))
hold on, plot(((Type(2).Data(5).t - Type(2).Data(5).t(end)) * 1000), numnansp2(1:length(Type(2).Data(5).t)), 'Color',linecolor(2,:))
set(gca,'Ylim', [0, 15], 'Xlim',[-100 0])


xlabel('Time (ms)', 'FontWeight', 'bold','FontSize',14)
ylabel('n', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca, 'box', 'off')
set(s2,'Position',[.13 .33 .7750 .1])
%% STEP 16: Box Plots for peak speed

linecolor = [[1 0 0] ;[0 0 1] ;[1 .5 0] ;[0 1 0]];


figure('Position',[300 300 500 500])
for p = 1:4
    index = [2 4 1 3];
    xdata = repmat(index(p), length(maxSpeed{1,p}), 1);
    scatter(xdata, maxSpeed{1,p}, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',linecolor(p,:),'MarkerEdgeColor',linecolor(p,:));
    hold on
    
    plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(maxSpeed{1,p}, 1), 2, 1), 'k-','LineWidth',3)
    hold on
    plot(repmat(xdata(1),2,1),[mean(maxSpeed{1,p}, 1)-(std(maxSpeed{1,p})/sqrt(length(maxSpeed{1,p})));mean(maxSpeed{1,p}, 1)+(std(maxSpeed{1,p})/sqrt(length(maxSpeed{1,p})))],'k-','LineWidth',1.5)
    hold on
end

ylabel('Peak Speed (m/s)', 'FontWeight', 'bold')
ylim([0 0.4])
xlim([.5 4.5])
ax = gca;
ax.FontSize = 14;

xticks([1 2 3 4])
yticks([0 .1 .2 .3 .4])
xticklabels({'GF1>+','GF1>Kir2.1','GF1>+','GF1>Kir2.1'})
xlabel('Capture                      Escape')
%% STEP 16B: Box plots for average speed

linecolor = [[1 0 0] ;[0 0 1] ;[1 .5 0] ;[0 1 0]];


figure('Position',[300 300 800 500])
for p = 1:4
    index = [2 4 1 3];
    xdata = repmat(index(p), length(avgSpeed{1,p}), 1);
    scatter(xdata, avgSpeed{1,p}, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',linecolor(p,:),'MarkerEdgeColor',linecolor(p,:));
    hold on
    
    plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(avgSpeed{1,p}, 1), 2, 1), 'k-','LineWidth',3)
    hold on
    plot(repmat(xdata(1),2,1),[mean(avgSpeed{1,p}, 1)-(std(avgSpeed{1,p})/sqrt(length(avgSpeed{1,p})));mean(avgSpeed{1,p}, 1)+(std(avgSpeed{1,p})/sqrt(length(avgSpeed{1,p})))],'k-','LineWidth',1.5)
    hold on
end

ylabel('Average Speed (m/s)', 'FontWeight', 'bold')
ylim([0 0.3])
xlim([.5 4.5])
ax = gca;
ax.FontSize = 14;

xticks([1 2 3 4])
yticks([0 .1 .2 .3])
xticklabels({'GF1>+','GF1>Kir2.1','GF1>+','GF1>Kir2.1'})
xlabel('Capture                                           Escape')
%% STEP 16C: Box plots for peak acceleration

linecolor = [[1 0 0] ;[0 0 1] ;[1 .5 0] ;[0 1 0]];

%DL Accel
figure('Position',[300 300 500 500])
for p = 1:4
    
    index = [2 4 1 3];
    xdata = repmat(index(p), length(peakAcc{1,p}), 1);
    scatter(xdata, peakAcc{1,p}, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',linecolor(p,:),'MarkerEdgeColor',linecolor(p,:));
    hold on
    
    plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(peakAcc{1,p}, 1), 2, 1), 'k-','LineWidth',3)
    hold on
    plot(repmat(xdata(1),2,1),[mean(peakAcc{1,p}, 1)-(std(peakAcc{1,p})/sqrt(length(peakAcc{1,p})));mean(peakAcc{1,p}, 1)+(std(peakAcc{1,p})/sqrt(length(peakAcc{1,p})))],'k-','LineWidth',1.5)
    hold on
end

%title('GF1 > +','FontSize',14)
ylabel('Peak Acceleration (m/s^{2})', 'FontWeight', 'bold')
ylim([0 6])
xlim([.5 4.5])
ax = gca;
ax.FontSize = 14;

xticks([1 2 3 4])

xticklabels({'GF1>+','GF1>Kir2.1','GF1>+','GF1>Kir2.1'})
xlabel('Capture                      Escape')

%% STEP 16D: Box plots for average acceleration

linecolor = [[1 0 0] ;[0 0 1] ;[1 .5 0] ;[0 1 0]];

%DL Accel
figure('Position',[300 300 800 500])
for p = 1:4
    
    index = [2 4 1 3];
    xdata = repmat(index(p), length(accAvg{1,p}), 1);
    scatter(xdata, accAvg{1,p}, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',linecolor(p,:),'MarkerEdgeColor',linecolor(p,:));
    hold on
    
    plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(accAvg{1,p}, 1), 2, 1), 'k-','LineWidth',3)
    hold on
    plot(repmat(xdata(1),2,1),[mean(accAvg{1,p}, 1)-(std(accAvg{1,p})/sqrt(length(accAvg{1,p})));mean(accAvg{1,p}, 1)+(std(accAvg{1,p})/sqrt(length(accAvg{1,p})))],'k-','LineWidth',1.5)
    hold on
end

%title('GF1 > +','FontSize',14)
ylabel('Average Acceleration (m/s^{2})', 'FontWeight', 'bold')
ylim([-3 3])
xlim([.5 4.5])
ax = gca;
ax.FontSize = 14;

xticks([1 2 3 4])

xticklabels({'GF1>+','GF1>Kir2.1','GF1>+','GF1>Kir2.1'})
xlabel('Capture                                           Escape')
%% STEP 17: Box plots for angular size
linecolor = [[1 0 0] ;[0 0 1] ;[1 .5 0] ;[0 1 0]];
thresh = 38.5;

%Average Acceleration
figure('Position',[300 300 500 500])
for p = 1:4
    
    index = [2 4 1 3];
    xdata = repmat(index(p), length(angSize{1,p}), 1);
    scatter(xdata, angSize{1,p}, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',linecolor(p,:),'MarkerEdgeColor',linecolor(p,:));
    hold on
    
    plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(angSize{1,p}, 1), 2, 1), 'k-','LineWidth',3)
    hold on
    plot(repmat(xdata(1),2,1),[mean(angSize{1,p}, 1)-(std(angSize{1,p})/sqrt(length(angSize{1,p})));mean(angSize{1,p}, 1)+(std(angSize{1,p})/sqrt(length(angSize{1,p})))],'k-','LineWidth',1.5)
    hold on
end

GFthresh = plot([0 5],[thresh thresh],'k--');


%title('GF1 > +','FontSize',12)
ylabel('Angular Size at Capture/Escape (ï¿½)', 'FontWeight', 'bold')
%xlabel('Fly Genotype', 'FontWeight', 'bold')
ylim([0 120])
xlim([.5 4.5])
ax = gca;
ax.FontSize = 14;
legend(GFthresh,{['Giant Fiber' 10 'Angular Threshold']},'Location',[0.73 0.9 0.1 0.1] )

% labels = {'line1 line2','line1 line2','line1 line2'};
% labels = cellfun(@(x) strrep(x,' ','\newline'), labels,'UniformOutput',false);
% a = gca;
% a.XTickLabel = labels;


xticks([1 2 3 4])
xticklabels({'GF1>+','GF1>Kir2.1','GF1>+','GF1>Kir2.1'})
xlabel('Capture                      Escape')
% set(gca,'Ydir','reverse')
% set(gca,'yscale','log')


% a2 = axes('YAxisLocation', 'Right');
% % Hide second plot.
% set(a2, 'color', 'none')
% set(a2, 'XTick', [])
% % Set scala for second Y.
% set(a2, 'YLim', [20 25])

%% STEP 18: distance from fly at TOC/Escape


linecolor = [[1 0 0] ;[0 0 1] ;[1 .5 0] ;[0 1 0]];

%Average Acceleration
figure('Position',[300 300 800 500])
for p = 1:4
    
    index = [2 4 1 3];
    xdata = repmat(index(p), length(distTOC{1,p}), 1);
    scatter(xdata, 1000*distTOC{1,p}, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',linecolor(p,:),'MarkerEdgeColor',linecolor(p,:));
    hold on
    
    plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], 1000*(repmat(mean(distTOC{1,p}, 1), 2, 1)), 'k-','LineWidth',3)
    hold on
    plot(repmat(xdata(1),2,1),1000*[mean(distTOC{1,p}, 1)-(std(distTOC{1,p})/sqrt(length(distTOC{1,p})));mean(distTOC{1,p}, 1)+(std(distTOC{1,p})/sqrt(length(distTOC{1,p})))],'k-','LineWidth',1.5)
    hold on
end


%title('GF1 > +','FontSize',12)
ylabel('Distance to Fly at Capture/Escape (mm)', 'FontWeight', 'bold','FontSize',14)
%ylim([y1min y1max])
ylim([0.5 100])
xlim([.5 4.5])
ax = gca;
ax.FontSize = 14;


xticks([1 2 3 4])
xticklabels({'GF1>+','GF1>Kir2.1','GF1>+','GF1>Kir2.1'})
xlabel('Capture                                           Escape')
% set(gca,'Ydir','reverse')
set(gca,'yscale','log')

%
% a2 = axes('YAxisLocation', 'Right');
% % Hide second plot.
% set(a2, 'color', 'none')
% set(a2, 'XTick', [])
% % Set scala for second Y.
% set(a2, 'YLim', [20 25])

%% STEP 19: Box plots for azimuith and elevation

linecolor = [[1 0 0] ;[0 0 1] ;[1 .5 0] ;[0 1 0]];

%Average Azimuth
figure('Position',[300 300 600 500])
for p = 1:4
    
    index = [2 4 1 3];
    xdata = repmat(index(p), length(avgAzimuth{1,p}), 1);
    scatter(xdata, avgAzimuth{1,p}, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',linecolor(p,:),'MarkerEdgeColor',linecolor(p,:));
    hold on
    
    plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(avgAzimuth{1,p}, 1), 2, 1), 'k-','LineWidth',3)
    hold on
    plot(repmat(xdata(1),2,1),[mean(avgAzimuth{1,p}, 1)-(std(avgAzimuth{1,p})/sqrt(length(avgAzimuth{1,p})));mean(avgAzimuth{1,p}, 1)+(std(avgAzimuth{1,p})/sqrt(length(avgAzimuth{1,p})))],'k-','LineWidth',1.5)
    hold on
end

%title('GF1 > +','FontSize',12)
ylabel('Average Azimuth Angle (ï¿½)', 'FontWeight', 'bold')
% xlabel('Fly Genotype', 'FontWeight', 'bold')
ylim([0 90])
xlim([.5 4.5])
ax = gca;
ax.FontSize = 14;

xticks([1 2 3 4])
xticklabels({'GF1>+','GF1>Kir2.1','GF1>+','GF1>Kir2.1'})
xlabel('Capture                             Escape')

% labels = {'line1 line2','line1 line2','line1 line2'};
% labels = cellfun(@(x) strrep(x,' ','\newline'), labels,'UniformOutput',false);
% a = gca;
% a.XTickLabel = labels;


linecolor = [[1 0 0] ;[0 0 1] ;[1 .5 0] ;[0 1 0]];

%Average Elevation
figure('Position',[300 300 600 500])
for p = 1:4
    
    index = [2 4 1 3];
    xdata = repmat(index(p), length(avgElevation{1,p}), 1);
    scatter(xdata, avgElevation{1,p}, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',linecolor(p,:),'MarkerEdgeColor',linecolor(p,:));
    hold on
    
    plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(avgElevation{1,p}, 1), 2, 1), 'k-','LineWidth',3)
    hold on
    plot(repmat(xdata(1),2,1),[mean(avgElevation{1,p}, 1)-(std(avgElevation{1,p})/sqrt(length(avgElevation{1,p})));mean(avgElevation{1,p}, 1)+(std(avgElevation{1,p})/sqrt(length(avgElevation{1,p})))],'k-','LineWidth',1.5)
    hold on
end

%title('GF1>+','FontSize',12)
ylabel('Average Elevation Angle (ï¿½)', 'FontWeight', 'bold')
% xlabel('Fly Genotype', 'FontWeight', 'bold')
ylim([-40 40])
xlim([.5 4.5])
ax = gca;
ax.FontSize = 14;

xticks([1 2 3 4])
xticklabels({'GF1>+','GF1>Kir2.1','GF1>+','GF1>Kir2.1'})
xlabel('Capture                             Escape')
%% STEP 20: Bar plots for velocity distribution
%minbin = 0.05;
%maxbin = 0.35;
%binsize = 0.15;
pctcapDL = [];
pctcapKir = [];

dlmaxspeed=[maxSpeed{3};maxSpeed{4}];
numbins = 4;
bins=min(dlmaxspeed):(max(dlmaxspeed)-min(dlmaxspeed))/numbins:max(dlmaxspeed);
minbin = bins(1);
maxbin = bins(end);
binsize = bins(2)-bins(1);

figure('Position',[300 300 700 500])
%for i = minbin:binsize:(maxbin-binsize)
for i = bins(1:end-1)
    for p = 1:4
        index= find(maxSpeed{1,p}>=i & maxSpeed{1,p}<(i+binsize));
        count(p) = length(index);
    end
    j = round(((i-minbin)/binsize)+1);
    pctcapDL(j) = count(3)/(count(3)+count(4));
    pctcapKir(j) = count(1)/(count(1)+count(2));
    countDL(j) = (count(3)+count(4));
    countKir(j) = (count(1)+count(2));
end
x = minbin+(binsize/2):binsize:maxbin-(binsize/2);
b = bar(x,100*[pctcapDL' pctcapKir'],'LineWidth',1.5);
b(1).FaceColor = [.5 .5 .5];
b(2).FaceColor = [1 1 1];

ylabel('Percent Capture (%)', 'FontWeight', 'bold')
xlabel('Peak Speed (m/s)', 'FontWeight', 'bold')
legend('GF1>+','GF1>Kir2.1','Location', 'northeast')
% ylim([0 0.6])
% xlim([2.5 4.5])
ax = gca;
ax.FontSize = 14;

%xticks(x)
set(gca, 'box', 'off')

figure('Position',[300 300 700 500])
x = minbin+(binsize/2):binsize:maxbin-(binsize/2);
b = plot(x,100*[pctcapDL' pctcapKir'],'LineWidth',1.5);
b(1).Color = [0 0 0];
b(2).Color = [.75 .75 .75];
b(1).Marker = 'o';
b(2).Marker = 'o';
b(1).MarkerFaceColor = [0 0 0];
b(2).MarkerFaceColor = [.75 .75 .75];

ylabel('Percent Capture (%)', 'FontWeight', 'bold')
xlabel('Peak Speed (m/s)', 'FontWeight', 'bold')
legend('GF1>+','GF1>Kir2.1','Location', 'northeast')
% ylim([0 0.6])
% xlim([2.5 4.5])
ax = gca;
ax.FontSize = 14;

%xticks(x)
set(gca, 'box', 'off')
%% STEP 21: Bar plots for acceleration distribution
% minbin = 7.5;
% maxbin = 32.5;
% binsize = 5;
minbin = 0;
maxbin = 6;
binsize = 2;

pctcapDL = [];
pctcapKir = [];
figure('Position',[300 300 700 500])
for i = minbin:binsize:(maxbin-binsize)
    for p = 1:4
        index= find(peakAcc{1,p}>=i & peakAcc{1,p}<(i+binsize));
        count(p) = length(index);
    end
    j = round(((i-minbin)/binsize)+1);
    pctcapDL(j) = count(4)/(count(3)+count(4));
    pctcapKir(j) = count(2)/(count(1)+count(2));
    countDL(j) = (count(3)+count(4));
    countKir(j) = (count(1)+count(2));
end
x = minbin+(binsize/2):binsize:maxbin-(binsize/2);
b = bar(x,100*[pctcapDL' pctcapKir'],'LineWidth',1.5);
b(1).FaceColor = [.5 .5 .5];
b(2).FaceColor = [1 1 1];

ylabel('Percent Escape (%)', 'FontWeight', 'bold')
xlabel('Peak Acceleration (m/s^{2})', 'FontWeight', 'bold')
legend('GF1>+','GF1>Kir2.1','Location', 'northeast')
% ylim([0 0.6])
% xlim([2.5 4.5])
ax = gca;
ax.FontSize = 14;

xticks(x)
set(gca, 'box', 'off')


figure('Position',[300 300 700 500])
for i = minbin:binsize:(maxbin-binsize)
    for p = 1:4
        index= find(peakAcc{1,p}>=i & peakAcc{1,p}<(i+binsize));
        count(p) = length(index);
    end
    j = round(((i-minbin)/binsize)+1);
    pctcapDL(j) = count(4)/(count(3)+count(4));
    pctcapKir(j) = count(2)/(count(1)+count(2));
    countDL(j) = (count(3)+count(4));
    countKir(j) = (count(1)+count(2));
end
x = minbin+(binsize/2):binsize:maxbin-(binsize/2);
b = plot(x,100*[pctcapDL' pctcapKir'],'LineWidth',1.5);
b(1).Color = [0 0 0];
b(2).Color = [.75 .75 .75];
b(1).Marker = 'o';
b(2).Marker = 'o';
b(1).MarkerFaceColor = [0 0 0];
b(2).MarkerFaceColor = [.75 .75 .75];

ylabel('Percent Escape (%)', 'FontWeight', 'bold')
xlabel('Peak Acceleration (m/s^{2})', 'FontWeight', 'bold')
legend('GF1>+','GF1>Kir2.1','Location', 'northeast')
% ylim([0 0.6])
% xlim([2.5 4.5])
ax = gca;
ax.FontSize = 14;

xticks(x)
set(gca, 'box', 'off')

%% STEP 20B: Line plots for velocity distribution - use mean
minbin = 0;
maxbin = .48;
binsize = .24;
pctcapDL = [];
pctcapKir = [];

figure('Position',[300 300 700 500])
for i = minbin:binsize:(maxbin-binsize)
    for p = 1:4
        index= find(maxSpeed{1,p}>=i & maxSpeed{1,p}<(i+binsize));
        count(p) = length(index);
    end
    j = round(((i-minbin)/binsize)+1);
    pctcapDL(j) = count(4)/(count(3)+count(4));
    pctcapKir(j) = count(2)/(count(1)+count(2));
    countDL(j) = (count(3)+count(4));
    countKir(j) = (count(1)+count(2));
end
x = minbin+(binsize/2):binsize:maxbin-(binsize/2);
b = plot(x,100*[pctcapDL' pctcapKir'],'LineWidth',1.5);
b(1).Color = [0 0 0];
b(2).Color = [.75 .75 .75];
b(1).Marker = 'o';
b(2).Marker = 'o';
b(1).MarkerFaceColor = [0 0 0];
b(2).MarkerFaceColor = [.75 .75 .75];

ylabel('Percent Escape (%)', 'FontWeight', 'bold')
xlabel('Peak Speed (m/s)', 'FontWeight', 'bold')
legend('GF1>+','GF1>Kir2.1','Location', 'northeast')
 ylim([0 100])
% xlim([2.5 4.5])
ax = gca;
ax.FontSize = 14;

xticks(x)
xticklabels({'<0.24','>0.24'})
set(gca, 'box', 'off')
%% STEP 21B: Line plots for acceleration distribution - use mean
minbin = 0;
maxbin = 3.8*2;
binsize = 3.8;

ylabel('Percent Escape (%)', 'FontWeight', 'bold')
xlabel('Peak Acceleration (m/s^{2})', 'FontWeight', 'bold')
legend('GF1>+','GF1>Kir2.1','Location', 'northeast')
ylim([0 100])
% xlim([2.5 4.5])
ax = gca;
ax.FontSize = 14;

xticks(x)
set(gca, 'box', 'off')


figure('Position',[300 300 700 500])
for i = minbin:binsize:(maxbin-binsize)
    for p = 1:4
        index= find(peakAcc{1,p}>=i & peakAcc{1,p}<(i+binsize));
        count(p) = length(index);
    end
    j = round(((i-minbin)/binsize)+1);
    pctcapDL(j) = count(4)/(count(3)+count(4));
    pctcapKir(j) = count(2)/(count(1)+count(2));
    countDL(j) = (count(3)+count(4));
    countKir(j) = (count(1)+count(2));
end
x = minbin+(binsize/2):binsize:maxbin-(binsize/2);
b = plot(x,100*[pctcapDL' pctcapKir'],'LineWidth',1.5);
b(1).Color = [0 0 0];
b(2).Color = [.75 .75 .75];
b(1).Marker = 'o';
b(2).Marker = 'o';
b(1).MarkerFaceColor = [0 0 0];
b(2).MarkerFaceColor = [.75 .75 .75];

ylabel('Percent Escape (%)', 'FontWeight', 'bold')
xlabel('Peak Acceleration (m/s^{2})', 'FontWeight', 'bold')
legend('GF1>+','GF1>Kir2.1','Location', 'northeast')
% ylim([0 0.6])
% xlim([2.5 4.5])
ax = gca;
ax.FontSize = 14;

xticks(x)
xticklabels({'<0.38','>0.38'})
set(gca, 'box', 'off')
%% STEP 22: Azimuth Plots Combined



linecolor = [[1 0 0] ;[0 0 1] ;[1 .5 0] ;[0 1 0]];
%Setup speed variables for plotting
speedp1 = padcat(flipud(Type(1).Data(1).damAzi), flipud(Type(1).Data(2).damAzi), flipud(Type(1).Data(3).damAzi), flipud(Type(1).Data(4).damAzi), flipud(Type(1).Data(5).damAzi)...
    ,flipud(Type(1).Data(6).damAzi), flipud(Type(1).Data(7).damAzi), flipud(Type(1).Data(8).damAzi), flipud(Type(1).Data(6).damAzi), flipud(Type(1).Data(10).damAzi), flipud(Type(1).Data(11).damAzi)...
    ,flipud(Type(1).Data(12).damAzi), flipud(Type(1).Data(13).damAzi));
speedp1 = flipud(speedp1); %concat speed matricies
numnansp1 = sum(~isnan(speedp1),2); %number of traces averaged as a function of time
meanp1Vec = nanmean(speedp1,2); %mean speed as a function of time
stdp1Vec = nanstd(speedp1, 0, 2); %standard deviation of speed

speedp2 = padcat(flipud(Type(2).Data(1).damAzi), flipud(Type(2).Data(2).damAzi), flipud(Type(2).Data(3).damAzi), flipud(Type(2).Data(4).damAzi), flipud(Type(2).Data(5).damAzi)...
    , flipud(Type(2).Data(6).damAzi));
speedp2 = flipud(speedp2);
numnansp2 = sum(~isnan(speedp2),2);
meanp2Vec = nanmean(speedp2,2);
stdp2Vec = nanstd(speedp2, 0, 2);

speedp3 = padcat(flipud(Type(3).Data(1).damAzi), flipud(Type(3).Data(2).damAzi), flipud(Type(3).Data(3).damAzi), flipud(Type(3).Data(4).damAzi), flipud(Type(3).Data(5).damAzi));
speedp3 = flipud(speedp3);
numnansp3 = sum(~isnan(speedp3),2);
meanp3Vec = nanmean(speedp3, 2);
stdp3Vec = nanstd(speedp3, 0, 2);

speedp4 = padcat(flipud(Type(4).Data(1).damAzi), flipud(Type(4).Data(2).damAzi), flipud(Type(4).Data(3).damAzi), flipud(Type(4).Data(4).damAzi), flipud(Type(4).Data(5).damAzi)...
    ,flipud(Type(4).Data(6).damAzi), flipud(Type(4).Data(7).damAzi));
speedp4 = flipud(speedp4);
numnansp4 = sum(~isnan(speedp4),2);
meanp4Vec = nanmean(speedp4,02);
stdp4Vec = nanstd(speedp4, 0, 2);

% Plot DL Speeds - Capture and Escape
figure('Position',[500 300 700 635])
s1 = subplot(2, 1, 1);
%
%                     plot((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000, meanp3Vec(1:length(Type(3).Data(1).t)),'Color',linecolor(3,:), 'LineWidth', 2); %mean
%                     plot((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000, (meanp3Vec(1:length(Type(3).Data(1).t)) - stdp3Vec(1:length(Type(3).Data(1).t))), '--','Color',linecolor(3,:)) %lower standard dev line
%                     plot((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000, (meanp3Vec(1:length(Type(3).Data(1).t)) + stdp3Vec(1:length(Type(3).Data(1).t))), '--','Color',linecolor(3,:)) %upper standard dev line
%
%                    plot((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000, meanp4Vec(1:length(Type(4).Data(6).t)), 'Color',linecolor(4,:), 'LineWidth', 2);
%                     plot((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000, (meanp4Vec(1:length(Type(4).Data(6).t)) - stdp4Vec(1:length(Type(4).Data(6).t))), '--','Color',linecolor(4,:))
%                     plot((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000, (meanp4Vec(1:length(Type(4).Data(6).t)) + stdp4Vec(1:length(Type(4).Data(6).t))), '--','Color',linecolor(4,:))
%
%                     plot((Type(1).Data(3).t - Type(1).Data(3).t(end))*1000, meanp1Vec(1:length(Type(1).Data(3).t)),'Color',linecolor(1,:), 'LineWidth', 2);
%                     plot((Type(1).Data(3).t - Type(1).Data(3).t(end))*1000, (meanp1Vec(1:length(Type(1).Data(3).t)) - stdp1Vec(1:length(Type(1).Data(3).t))), '--','Color',linecolor(1,:))
%                     plot((Type(1).Data(3).t - Type(1).Data(3).t(end))*1000, (meanp1Vec(1:length(Type(1).Data(3).t)) + stdp1Vec(1:length(Type(1).Data(3).t))), '--r','Color',linecolor(1,:))
%
%                    plot((Type(2).Data(5).t - Type(2).Data(5).t(end))*1000, meanp2Vec(1:length(Type(2).Data(5).t)),'Color',linecolor(2,:), 'LineWidth', 2);
%                     plot((Type(2).Data(5).t - Type(2).Data(5).t(end))*1000, (meanp2Vec(1:length(Type(2).Data(5).t)) - stdp2Vec(1:length(Type(2).Data(5).t))), '--','Color',linecolor(2,:))
%                     plot((Type(2).Data(5).t - Type(2).Data(5).t(end))*1000, (meanp2Vec(1:length(Type(2).Data(5).t)) + stdp2Vec(1:length(Type(2).Data(5).t))), '--b','Color',linecolor(2,:))
%
three = shadedErrorBar(((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000), rad2deg(meanp3Vec(1:length(Type(3).Data(1).t))), rad2deg(stdp3Vec(1:length(Type(3).Data(1).t))),'lineprops',{'Color',linecolor(3,:),'LineWidth',2},'transparent', 1); %shaded area for standard dev
hold on
four = shadedErrorBar(((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000), rad2deg(meanp4Vec(1:length(Type(4).Data(6).t))), rad2deg(stdp4Vec(1:length(Type(4).Data(6).t))), 'lineprops',{'Color', linecolor(4,:),'LineWidth',2}, 'transparent', 1);
hold on
one = shadedErrorBar((Type(1).Data(3).t - Type(1).Data(3).t(end))*1000, rad2deg(meanp1Vec(1:length(Type(1).Data(3).t))), rad2deg(stdp1Vec(1:length(Type(1).Data(3).t))), 'lineprops',{'Color', linecolor(1,:),'LineWidth',2}, 'transparent', 1);
hold on
two = shadedErrorBar((Type(2).Data(5).t - Type(2).Data(5).t(end))*1000, rad2deg(meanp2Vec(1:length(Type(2).Data(5).t))), rad2deg(stdp2Vec(1:length(Type(2).Data(5).t))), 'lineprops',{'Color', linecolor(2,:),'LineWidth',2}, 'transparent', 1);


%xlabel('Time (ms)', 'FontWeight', 'bold','FontSize',12)
ylabel('Azimuth Angle (ï¿½)', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[0,120],'Xlim',[-100 0])
set(s1,'Position',[.13 .5 .7750 .45])
legend([one.mainLine,two.mainLine,three.mainLine,four.mainLine], 'GF1>Kir2.1 - Capture', 'GF1>Kir2.1 - Escape', 'GF1>+ - Capture', 'GF1>+ - Escape', 'Location', 'northwest', 'Orientation', 'horizontal');
set(legend, 'NumColumns' ,2)
set(gca, 'box', 'off')
hold off

s2 = subplot(2, 1, 2);
plot(((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000), numnansp3(1:length(Type(3).Data(1).t)), 'Color',linecolor(3,:)) %number of traces averaged
hold on, plot(((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000), numnansp4(1:length(Type(4).Data(6).t)), 'Color',linecolor(4,:))
hold on, plot(((Type(1).Data(3).t - Type(1).Data(3).t(end)) * 1000), numnansp1(1:length(Type(1).Data(3).t)), 'Color',linecolor(1,:))
hold on, plot(((Type(2).Data(5).t - Type(2).Data(5).t(end)) * 1000), numnansp2(1:length(Type(2).Data(5).t)), 'Color',linecolor(2,:))
set(gca,'Ylim', [0, 15], 'Xlim',[-100 0])


xlabel('Time (ms)', 'FontWeight', 'bold','FontSize',14)
ylabel('n', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca, 'box', 'off')
set(s2,'Position',[.13 .33 .7750 .1])

%% STEP 30:Individual elevation

%Setup speed variables for plotting
speedp1 = padcat(flipud(Type(1).Data(1).speedVec_Dams), flipud(Type(1).Data(2).speedVec_Dams), flipud(Type(1).Data(3).speedVec_Dams), flipud(Type(1).Data(4).speedVec_Dams), flipud(Type(1).Data(5).speedVec_Dams)...
    ,flipud(Type(1).Data(6).speedVec_Dams), flipud(Type(1).Data(7).speedVec_Dams), flipud(Type(1).Data(8).speedVec_Dams), flipud(Type(1).Data(6).speedVec_Dams), flipud(Type(1).Data(10).speedVec_Dams), flipud(Type(1).Data(11).speedVec_Dams)...
    ,flipud(Type(1).Data(12).speedVec_Dams), flipud(Type(1).Data(13).speedVec_Dams));
speedp1 = flipud(speedp1); %concat speed matricies
numnansp1 = sum(~isnan(speedp1),2); %number of traces averaged as a function of time
meanp1Vec = nanmean(speedp1, 2); %mean speed as a function of time
stdp1Vec = nanstd(speedp1, 0, 2); %standard deviation of speed

indlow = find(maxSpeed{1}<mean(maxSpeed{1}));
indhig = find(maxSpeed{1}>mean(maxSpeed{1}));

speedl = speedp1(:,indlow);
meanlVec = nanmean(speedl,2);
stdlVec = nanstd(speedl,0,2);

speedh = speedp1(:,indhig);
meanhVec = nanmean(speedh,2);
stdhVec = nanstd(speedh,0,2);

speedp2 = padcat(flipud(Type(2).Data(1).speedVec_Dams), flipud(Type(2).Data(2).speedVec_Dams), flipud(Type(2).Data(3).speedVec_Dams), flipud(Type(2).Data(4).speedVec_Dams), flipud(Type(2).Data(5).speedVec_Dams)...
    , flipud(Type(2).Data(6).speedVec_Dams));
speedp2 = flipud(speedp2);
numnansp2 = sum(~isnan(speedp2),2);
meanp2Vec = nanmean(speedp2, 2);
stdp2Vec = nanstd(speedp2, 0, 2);

speedp3 = padcat(flipud(Type(3).Data(1).speedVec_Dams), flipud(Type(3).Data(2).speedVec_Dams), flipud(Type(3).Data(3).speedVec_Dams), flipud(Type(3).Data(4).speedVec_Dams), flipud(Type(3).Data(5).speedVec_Dams));
speedp3 = flipud(speedp3);
numnansp3 = sum(~isnan(speedp3),2);
meanp3Vec = nanmean(speedp3, 2);
stdp3Vec = nanstd(speedp3, 0, 2);

speedp4 = padcat(flipud(Type(4).Data(1).speedVec_Dams), flipud(Type(4).Data(2).speedVec_Dams), flipud(Type(4).Data(3).speedVec_Dams), flipud(Type(4).Data(4).speedVec_Dams), flipud(Type(4).Data(5).speedVec_Dams)...
    ,flipud(Type(4).Data(6).speedVec_Dams), flipud(Type(4).Data(7).speedVec_Dams));
speedp4 = flipud(speedp4);
numnansp4 = sum(~isnan(speedp4),2);
meanp4Vec = nanmean(speedp4, 2);
stdp4Vec = nanstd(speedp4, 0, 2);

linecolor = [[1 0 0] ;[0 0 1] ;[1 .5 0] ;[0 1 0]];

% Plot DL - Capture Speed
figure('Position',[500 300 500 500])

% Uncomment for raw traces
for d = 1:length(Type(3).Data)
    
    plot((Type(3).Data(d).t - Type(3).Data(d).t(end)) * 1000,Type(3).Data(d).speedVec_Dams(1:length(Type(3).Data(d).t)), 'Color',linecolor(3,:));
    hold on
end

shadedErrorBar(((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000), meanp3Vec(1:length(Type(3).Data(1).t)), stdp3Vec(1:length(Type(3).Data(1).t)),'lineprops',{'Color',linecolor(3,:),'LineWidth',2},'transparent', 1); %shaded area for standard dev

for d = 1:length(Type(3).Data)
    time = (Type(3).Data(d).t - Type(3).Data(d).t(end)) * 1000;
    speed = Type(3).Data(d).speedVec_Dams(1:length(Type(3).Data(d).t));
    [m,ind] = max(speed);
    maxtime(d) = time(ind);
    plot((Type(3).Data(d).t - Type(3).Data(d).t(end)) * 1000,Type(3).Data(d).speedVec_Dams(1:length(Type(3).Data(d).t)), 'Color',linecolor(3,:));
    hold on
end


%title('GF1>+ - Capture','FontSize',14)
ylabel('Speed (m/s)', 'FontWeight', 'bold','FontSize',14)
xlabel('time(ms)', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[0,0.325],'Xlim',[-100 0])
hold on 
plot([-29,-29],[0,.3],'-k')
yticks([0 .1 .2 .3 .4])
set(gca, 'box', 'off')
%ax1 = gca;                   % gca = get current axis
%ax1.YAxis.Visible = 'off';   % remove y-axis
%ax1.XAxis.Visible = 'off';
hold off

% Plot Kir - Capture Speed
figure('Position',[500 300 500 500])

% Uncomment for raw traces
for d = 1:length(Type(1).Data)
    
    plot((Type(1).Data(d).t - Type(1).Data(d).t(end)) * 1000,Type(1).Data(d).speedVec_Dams(1:length(Type(1).Data(d).t)), 'Color',linecolor(1,:));
    hold on
end

shadedErrorBar(((Type(1).Data(3).t - Type(1).Data(3).t(end)) * 1000), meanp1Vec(1:length(Type(1).Data(3).t)), stdp1Vec(1:length(Type(1).Data(3).t)),'lineprops',{'Color',linecolor(1,:),'LineWidth',2},'transparent', 1); %shaded area for standard dev


%title('GF1>Kir2.1 - Capture','FontSize',14)
ylabel('Speed (m/s)', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[0,0.325],'Xlim',[-100 0])
set(gca, 'box', 'off')
ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';
yticks([0 .1 .2 .3 .4])

% Plot DL - Escape Speed
figure('Position',[500 300 500 500])

% Uncomment for raw traces
for d = 1:length(Type(4).Data)
    
    plot((Type(4).Data(d).t - Type(4).Data(d).t(end)) * 1000,Type(4).Data(d).speedVec_Dams(1:length(Type(4).Data(d).t)), 'Color',linecolor(4,:));
    hold on
end

shadedErrorBar(((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000), meanp4Vec(1:length(Type(4).Data(6).t)), stdp4Vec(1:length(Type(4).Data(6).t)),'lineprops',{'Color',linecolor(4,:),'LineWidth',2},'transparent', 1); %shaded area for standard dev


%title('GF1>+ - Escape','FontSize',14)
ylabel('Speed (m/s)', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[0,0.325],'Xlim',[-100 0])
set(gca, 'box', 'off')
ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';
yticks([0 .1 .2 .3 .4])
hold off

% Plot Kir - Escape Speed
figure('Position',[500 300 500 500])

% Uncomment for raw traces
for d = 1:length(Type(2).Data)
    
    plot((Type(2).Data(d).t - Type(2).Data(d).t(end)) * 1000,Type(2).Data(d).speedVec_Dams(1:length(Type(2).Data(d).t)), 'Color',linecolor(2,:));
    hold on
end

shadedErrorBar(((Type(2).Data(5).t - Type(2).Data(5).t(end)) * 1000), meanp2Vec(1:length(Type(2).Data(5).t)), stdp2Vec(1:length(Type(2).Data(5).t)),'lineprops',{'Color',linecolor(2,:),'LineWidth',2},'transparent', 1); %shaded area for standard dev


%title('GF1>Kir2.1 - Escape','FontSize',14)
ylabel('Speed (m/s)', 'FontWeight', 'bold','FontSize',14)
xlabel('Time (ms)', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[0,0.325],'Xlim',[-100 0])
set(gca, 'box', 'off')
ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';
yticks([0 .1 .2 .3 .4])

%% STEP 23: Elevation Plots Combined

linecolor = [[1 0 0] ;[0 0 1] ;[1 .5 0] ;[0 1 0]];
%Setup speed variables for plotting
speedp1 = padcat(flipud(Type(1).Data(1).damEle), flipud(Type(1).Data(2).damEle), flipud(Type(1).Data(3).damEle), flipud(Type(1).Data(4).damEle), flipud(Type(1).Data(5).damEle)...
    ,flipud(Type(1).Data(6).damEle), flipud(Type(1).Data(7).damEle), flipud(Type(1).Data(8).damEle), flipud(Type(1).Data(6).damEle), flipud(Type(1).Data(10).damEle), flipud(Type(1).Data(11).damEle)...
    ,flipud(Type(1).Data(12).damEle), flipud(Type(1).Data(13).damEle));
speedp1 = flipud(speedp1); %concat speed matricies
numnansp1 = sum(~isnan(speedp1),2); %number of traces averaged as a function of time
meanp1Vec = nanmean(speedp1,2); %mean speed as a function of time
stdp1Vec = nanstd(speedp1, 0, 2); %standard deviation of speed

speedp2 = padcat(flipud(Type(2).Data(1).damEle), flipud(Type(2).Data(2).damEle), flipud(Type(2).Data(3).damEle), flipud(Type(2).Data(4).damEle), flipud(Type(2).Data(5).damEle)...
    , flipud(Type(2).Data(6).damEle));
speedp2 = flipud(speedp2);
numnansp2 = sum(~isnan(speedp2),2);
meanp2Vec = nanmean(speedp2,2);
stdp2Vec = nanstd(speedp2, 0, 2);

speedp3 = padcat(flipud(Type(3).Data(1).damEle), flipud(Type(3).Data(2).damEle), flipud(Type(3).Data(3).damEle), flipud(Type(3).Data(4).damEle), flipud(Type(3).Data(5).damEle));
speedp3 = flipud(speedp3);
numnansp3 = sum(~isnan(speedp3),2);
meanp3Vec = nanmean(speedp3, 2);
stdp3Vec = nanstd(speedp3, 0, 2);

speedp4 = padcat(flipud(Type(4).Data(1).damEle), flipud(Type(4).Data(2).damEle), flipud(Type(4).Data(3).damEle), flipud(Type(4).Data(4).damEle), flipud(Type(4).Data(5).damEle)...
    ,flipud(Type(4).Data(6).damEle), flipud(Type(4).Data(7).damEle));
speedp4 = flipud(speedp4);
numnansp4 = sum(~isnan(speedp4),2);
meanp4Vec = nanmean(speedp4,02);
stdp4Vec = nanstd(speedp4, 0, 2);

% Plot DL Speeds - Capture and Escape
figure('Position',[500 300 700 635])
s1 = subplot(2, 1, 1);
%
%                     plot((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000, meanp3Vec(1:length(Type(3).Data(1).t)),'Color',linecolor(3,:), 'LineWidth', 2); %mean
%                     plot((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000, (meanp3Vec(1:length(Type(3).Data(1).t)) - stdp3Vec(1:length(Type(3).Data(1).t))), '--','Color',linecolor(3,:)) %lower standard dev line
%                     plot((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000, (meanp3Vec(1:length(Type(3).Data(1).t)) + stdp3Vec(1:length(Type(3).Data(1).t))), '--','Color',linecolor(3,:)) %upper standard dev line
%
%                    plot((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000, meanp4Vec(1:length(Type(4).Data(6).t)), 'Color',linecolor(4,:), 'LineWidth', 2);
%                     plot((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000, (meanp4Vec(1:length(Type(4).Data(6).t)) - stdp4Vec(1:length(Type(4).Data(6).t))), '--','Color',linecolor(4,:))
%                     plot((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000, (meanp4Vec(1:length(Type(4).Data(6).t)) + stdp4Vec(1:length(Type(4).Data(6).t))), '--','Color',linecolor(4,:))
%
%                     plot((Type(1).Data(3).t - Type(1).Data(3).t(end))*1000, meanp1Vec(1:length(Type(1).Data(3).t)),'Color',linecolor(1,:), 'LineWidth', 2);
%                     plot((Type(1).Data(3).t - Type(1).Data(3).t(end))*1000, (meanp1Vec(1:length(Type(1).Data(3).t)) - stdp1Vec(1:length(Type(1).Data(3).t))), '--','Color',linecolor(1,:))
%                     plot((Type(1).Data(3).t - Type(1).Data(3).t(end))*1000, (meanp1Vec(1:length(Type(1).Data(3).t)) + stdp1Vec(1:length(Type(1).Data(3).t))), '--r','Color',linecolor(1,:))
%
%                    plot((Type(2).Data(5).t - Type(2).Data(5).t(end))*1000, meanp2Vec(1:length(Type(2).Data(5).t)),'Color',linecolor(2,:), 'LineWidth', 2);
%                     plot((Type(2).Data(5).t - Type(2).Data(5).t(end))*1000, (meanp2Vec(1:length(Type(2).Data(5).t)) - stdp2Vec(1:length(Type(2).Data(5).t))), '--','Color',linecolor(2,:))
%                     plot((Type(2).Data(5).t - Type(2).Data(5).t(end))*1000, (meanp2Vec(1:length(Type(2).Data(5).t)) + stdp2Vec(1:length(Type(2).Data(5).t))), '--b','Color',linecolor(2,:))
%
three = shadedErrorBar(((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000), rad2deg(meanp3Vec(1:length(Type(3).Data(1).t))), rad2deg(stdp3Vec(1:length(Type(3).Data(1).t))),'lineprops',{'Color',linecolor(3,:),'LineWidth',2},'transparent', 1); %shaded area for standard dev
hold on
four = shadedErrorBar(((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000), rad2deg(meanp4Vec(1:length(Type(4).Data(6).t))), rad2deg(stdp4Vec(1:length(Type(4).Data(6).t))), 'lineprops',{'Color', linecolor(4,:),'LineWidth',2}, 'transparent', 1);
hold on
one = shadedErrorBar((Type(1).Data(3).t - Type(1).Data(3).t(end))*1000, rad2deg(meanp1Vec(1:length(Type(1).Data(3).t))), rad2deg(stdp1Vec(1:length(Type(1).Data(3).t))), 'lineprops',{'Color', linecolor(1,:),'LineWidth',2}, 'transparent', 1);
hold on
two = shadedErrorBar((Type(2).Data(5).t - Type(2).Data(5).t(end))*1000, rad2deg(meanp2Vec(1:length(Type(2).Data(5).t))), rad2deg(stdp2Vec(1:length(Type(2).Data(5).t))), 'lineprops',{'Color', linecolor(2,:),'LineWidth',2}, 'transparent', 1);


%xlabel('Time (ms)', 'FontWeight', 'bold','FontSize',12)
ylabel('Elevation Angle (ï¿½)', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[-90,90],'Xlim',[-100 0])
set(s1,'Position',[.13 .5 .7750 .45])
legend([one.mainLine,two.mainLine,three.mainLine,four.mainLine], 'GF1>Kir2.1 - Capture', 'GF1>Kir2.1 - Escape', 'GF1>+ - Capture', 'GF1>+ - Escape', 'Location', 'northwest', 'Orientation', 'horizontal');
set(legend, 'NumColumns' ,2)
set(gca, 'box', 'off')
hold off

s2 = subplot(2, 1, 2);
plot(((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000), numnansp3(1:length(Type(3).Data(1).t)), 'Color',linecolor(3,:)) %number of traces averaged
hold on, plot(((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000), numnansp4(1:length(Type(4).Data(6).t)), 'Color',linecolor(4,:))
hold on, plot(((Type(1).Data(3).t - Type(1).Data(3).t(end)) * 1000), numnansp1(1:length(Type(1).Data(3).t)), 'Color',linecolor(1,:))
hold on, plot(((Type(2).Data(5).t - Type(2).Data(5).t(end)) * 1000), numnansp2(1:length(Type(2).Data(5).t)), 'Color',linecolor(2,:))
set(gca,'Ylim', [0, 15], 'Xlim',[-100 0])


xlabel('Time (ms)', 'FontWeight', 'bold','FontSize',14)
ylabel('n', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca, 'box', 'off')
set(s2,'Position',[.13 .33 .7750 .1])

%% STEP 24: Headsize Plot
headsizes = [3.1994 2.8572 3.1649 3.1661 3.2425 2.9415 2.9619 2.9619 3.1862 3.0545 3.017 3.035 2.8935 3.0679 3.2927 2.8469 3.0549 2.7872 2.9892 2.7194 2.8938 2.8996 2.6701 2.6069 2.7273 2.8049 3.0179 3.0187 2.8909 2.8565 2.9359 2.6738 2.9673];

%DL Accel
figure('Position',[300 300 300 500])


xdata = repmat(1, length(headsizes), 1);
scatter(xdata, headsizes, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor',[.75 .75 .75]);
hold on

plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(headsizes), 1, 2, 1), 'k-','LineWidth',3)
hold on
plot(repmat(xdata(1),2,1),[mean(headsizes, 2)-(std(headsizes)/sqrt(length(headsizes)));mean(headsizes, 2)+(std(headsizes)/sqrt(length(headsizes)))],'k-','LineWidth',1.5)
hold on

%title('GF1 > +','FontSize',12)
ylabel('Damselfly Head Width (mm)', 'FontWeight', 'bold')
ylim([2.5 3.5])
xlim([.5 1.5])
ax = gca;
ax.FontSize = 14;

xticks([5])

xticklabels({''})

%% STEP 25: Left-Right Plot

wt1 = [1 0.14 -0.17 -0.33 -0.29 -0.11 0.05 -0.6 -0.47];
kir1 = [0.13 0 0.12 0.12 0 -0.33];
kir1gal80 = [0.30 -0.05 0.04 0.06 -0.07 -0.03 0.38];
wt2 = [0.1 0.39 0.09 -0.03];
kir2 = [0.23 0.04 -0.56 -0.08];


figure('Position',[300 300 600 550])

xdata = repmat(1,length(wt1),1);
scatter(xdata, wt1, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor',[.75 .75 .75]);
hold on
plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(wt1), 2, 1), 'k-','LineWidth',3)
hold on
plot(repmat(xdata(1),2,1),[mean(wt1, 2)-(std(wt1)/sqrt(length(wt1)));mean(wt1, 2)+(std(wt1)/sqrt(length(wt1)))],'k-','LineWidth',1.5)
hold on

xdata = repmat(2,length(kir1),1);
scatter(xdata, kir1, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor',[.75 .75 .75]);
hold on
plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(kir1), 2, 1), 'k-','LineWidth',3)
hold on
plot(repmat(xdata(1),2,1),[mean(kir1, 2)-(std(kir1)/sqrt(length(kir1)));mean(kir1, 2)+(std(kir1)/sqrt(length(kir1)))],'k-','LineWidth',1.5)
hold on

xdata = repmat(3,length(kir1gal80),1);
scatter(xdata, kir1gal80, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor',[.75 .75 .75]);
hold on
plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(kir1gal80), 2, 1), 'k-','LineWidth',3)
hold on
plot(repmat(xdata(1),2,1),[mean(kir1gal80, 2)-(std(kir1gal80)/sqrt(length(kir1gal80)));mean(kir1gal80, 2)+(std(kir1gal80)/sqrt(length(kir1gal80)))],'k-','LineWidth',1.5)
hold on

xdata = repmat(4,length(wt2),1);
scatter(xdata,  wt2, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor',[.75 .75 .75]);
hold on
plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(wt2), 2, 1), 'k-','LineWidth',3)
hold on
plot(repmat(xdata(1),2,1),[mean(wt2, 2)-(std(wt2)/sqrt(length(wt2)));mean(wt2, 2)+(std(wt2)/sqrt(length(wt2)))],'k-','LineWidth',1.5)
hold on

xdata = repmat(5,length(kir2),1);
scatter(xdata, kir2, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor',[.75 .75 .75]);
hold on
plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(kir2), 2, 1), 'k-','LineWidth',3)
hold on
plot(repmat(xdata(1),2,1),[mean(kir2, 2)-(std(kir2)/sqrt(length(kir2)));mean(kir2, 2)+(std(kir2)/sqrt(length(kir2)))],'k-','LineWidth',1.5)
hold on

plot([0 6],[0 0],'k--')

ylabel('Wing Clipping Bias Index', 'FontWeight', 'bold')
ylim([-.8 1])
xlim([.5 5.5])
ax = gca;
ax.FontSize = 14;

xticks([1 2 3 4 5])
xticklabels({'GF1>+','GF1>Kir2.1','GF1>Kir2.1;Gal80ts','GF2>+','GF2>Kir2.1'})
xtickangle(30)




%% STEP 26: Prey Consumption Index Plots

pooled = [wt1 kir1 kir1gal80 wt2 kir2];
match2 = [0.31 0.29 0.07 0.50 0.68 0.39 0.75 0.45 1.00 0.13 0.33 0.45 -0.08 0.25 0.43];
match3 = [0.23 0.29 0.71 0.28 0.38 0.33 0.75 1.00];
match1 = [0.19 0.50 0.40 0.04 0.30 0.16 0.04 0.60];


figure('Position',[300 300 250 515])

xdata = repmat(1,length(pooled),1);
scatter(xdata, pooled, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor',[.75 .75 .75]);
hold on
plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(pooled), 2, 1), 'k-','LineWidth',3)
hold on
plot(repmat(xdata(1),2,1),[mean(pooled, 2)-(std(wt1)/sqrt(length(pooled)));mean(pooled, 2)+(std(wt1)/sqrt(length(pooled)))],'k-','LineWidth',1.5)
hold on

plot([0 2],[0 0],'k--')

%xlabel('A', 'FontWeight', 'bold')
ylabel('Wing Clip Bias Index', 'FontWeight', 'bold')
ylim([-.7 1])
xlim([.5 1.5])
ax = gca;
ax.FontSize = 14;

xticks(1)

labels = {' Pooled\Controls'};
labels = cellfun(@(x) strrep(x,'\','\newline'), labels,'UniformOutput',false);
% a = gca;
% a.XTickLabel = labels;

xticklabels(labels)
xtickangle(30)
ax.XRuler.TickLabelGapOffset = 0;

figure('Position',[300 300 600 600])

xdata = repmat(1,length(match1),1);
scatter(xdata, match1, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor',[.75 .75 .75]);
hold on
plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(match1), 2, 1), 'k-','LineWidth',3)
hold on
plot(repmat(xdata(1),2,1),[mean(match1, 2)-(std(match1)/sqrt(length(match1)));mean(match1, 2)+(std(match1)/sqrt(length(match1)))],'k-','LineWidth',1.5)
hold on

xdata = repmat(2,length(match2),1);
scatter(xdata, match2, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor',[.75 .75 .75]);
hold on
plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(match2), 2, 1), 'k-','LineWidth',3)
hold on
plot(repmat(xdata(1),2,1),[mean(match2, 2)-(std(match2)/sqrt(length(match2)));mean(match2, 2)+(std(match2)/sqrt(length(match2)))],'k-','LineWidth',1.5)
hold on

xdata = repmat(3,length(match3),1);
scatter(xdata,  match3, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor',[.75 .75 .75]);
hold on
plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(match3), 2, 1), 'k-','LineWidth',3)
hold on
plot(repmat(xdata(1),2,1),[mean(match3, 2)-(std(match3)/sqrt(length(match3)));mean(match3, 2)+(std(match3)/sqrt(length(match3)))],'k-','LineWidth',1.5)
hold on



plot([0 4],[0 0],'k--')

xlabel('Fly Genotype Matchup', 'FontWeight', 'bold','HorizontalAlignment','right')
ylabel('Prey Consumption Index', 'FontWeight', 'bold')
ylim([-.7 1])
xlim([.5 3.5])
ax = gca;
ax.FontSize = 14;

xticks([1 2 3])

labels = {'                      GF1>+\  GF1>Kir2.1;Gal80ts','         GF1>+\  GF1>Kir2.1','        GF2>+\  GF2>Kir2.1'};
labels = cellfun(@(x) strrep(x,'\','\newline'), labels,'UniformOutput',false);
% a = gca;
% a.XTickLabel = labels;

xticklabels(labels)
xtickangle(30)
ax.XRuler.TickLabelGapOffset = 0;

%% STEP 27: Num Eaten Plots

%GF1>+

left =  [9 8 5 4 10 8 10 2 4];
right = [0 6 7 8 18 10 9 8 11];

figure('Position',[300 300 200 300])
xdata = repmat(1,length(left),1);
scatter(xdata,  left, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor',[.75 .75 .75]);
hold on
plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(left), 2, 1), 'k-','LineWidth',3)
hold on
plot(repmat(xdata(1),2,1),[mean(left, 2)-(std(left)/sqrt(length(left)));mean(left, 2)+(std(left)/sqrt(length(left)))],'k-','LineWidth',1.5)
hold on

xdata = repmat(2,length(right),1);
scatter(xdata,  right, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor',[.75 .75 .75]);
hold on
plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(right), 2, 1), 'k-','LineWidth',3)
hold on
plot(repmat(xdata(1),2,1),[mean(right, 2)-(std(right)/sqrt(length(right)));mean(right, 2)+(std(right)/sqrt(length(right)))],'k-','LineWidth',1.5)
hold on

xticks([1 2])
ylabel('Number Eaten', 'FontWeight', 'bold')
ylim([0 20])
xlim([0.5 2.5])
xticklabels(["",""])


%GF1>Kir

left =  [9 10 14 14 10 3];
right = [7 10 11 11 10 6];

figure('Position',[300 300 200 300])
xdata = repmat(1,length(left),1);
scatter(xdata,  left, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor',[.75 .75 .75]);
hold on
plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(left), 2, 1), 'k-','LineWidth',3)
hold on
plot(repmat(xdata(1),2,1),[mean(left, 2)-(std(left)/sqrt(length(left)));mean(left, 2)+(std(left)/sqrt(length(left)))],'k-','LineWidth',1.5)
hold on

xdata = repmat(2,length(right),1);
scatter(xdata,  right, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor',[.75 .75 .75]);
hold on
plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(right), 2, 1), 'k-','LineWidth',3)
hold on
plot(repmat(xdata(1),2,1),[mean(right, 2)-(std(right)/sqrt(length(right)));mean(right, 2)+(std(right)/sqrt(length(right)))],'k-','LineWidth',1.5)
hold on

xticks([1 2])
%ylabel('Number Eaten', 'FontWeight', 'bold')
ylim([0 20])
xlim([0.5 2.5])
xticklabels(["",""])
yticklabels(["","","","",""])


%L=GF1>+, R=GF1>Kir

left =  [7  2  0  7  7 6];
right = [16 14 12 14 6 15];

figure('Position',[300 300 200 300])
xdata = repmat(1,length(left),1);
scatter(xdata,  left, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor',[.75 .75 .75]);
hold on
plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(left), 2, 1), 'k-','LineWidth',3)
hold on
plot(repmat(xdata(1),2,1),[mean(left, 2)-(std(left)/sqrt(length(left)));mean(left, 2)+(std(left)/sqrt(length(left)))],'k-','LineWidth',1.5)
hold on

xdata = repmat(2,length(right),1);
scatter(xdata,  right, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor',[.75 .75 .75]);
hold on
plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(right), 2, 1), 'k-','LineWidth',3)
hold on
plot(repmat(xdata(1),2,1),[mean(right, 2)-(std(right)/sqrt(length(right)));mean(right, 2)+(std(right)/sqrt(length(right)))],'k-','LineWidth',1.5)
hold on

xticks([1 2])
%ylabel('Number Eaten', 'FontWeight', 'bold')
ylim([0 20])
xlim([0.5 2.5])
xticklabels(["",""])

yticklabels(["","","","",""])

%L=GF1>Kir, R=GF1>+

left =  [16 8 13 8 15];
right = [3 3 10 3 9];

figure('Position',[300 300 200 300])
xdata = repmat(1,length(left),1);
scatter(xdata,  left, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor',[.75 .75 .75]);
hold on
plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(left), 2, 1), 'k-','LineWidth',3)
hold on
plot(repmat(xdata(1),2,1),[mean(left, 2)-(std(left)/sqrt(length(left)));mean(left, 2)+(std(left)/sqrt(length(left)))],'k-','LineWidth',1.5)
hold on

xdata = repmat(2,length(right),1);
scatter(xdata,  right, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor',[.75 .75 .75]);
hold on
plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(right), 2, 1), 'k-','LineWidth',3)
hold on
plot(repmat(xdata(1),2,1),[mean(right, 2)-(std(right)/sqrt(length(right)));mean(right, 2)+(std(right)/sqrt(length(right)))],'k-','LineWidth',1.5)
hold on

xticks([1 2])
%ylabel('Number Eaten', 'FontWeight', 'bold')
ylim([0 20])
xlim([0.5 2.5])
xticklabels(["",""])
yticklabels(["","","","",""])

%% STEP 28: Pez TO and Short

% Order is xDL, xWeak Kir, xStrong Kir


load(fullfile(pwd, 'FlyPEZ Data\0038000015110253_dataForVisualization.mat'))
load(fullfile(pwd, 'FlyPEZ Data\0038000015110253_manualAnnotations.mat'))
graphTable1 = join(graphTable,manualAnnotations,'Keys','Row');
load(fullfile(pwd, 'FlyPEZ Data\0038000015110321_dataForVisualization.mat')) 
load(fullfile(pwd, 'FlyPEZ Data\0038000015110321_manualAnnotations.mat'))
graphTable2 = join(graphTable,manualAnnotations,'Keys','Row');
GF1DL10 = vertcat(graphTable1,graphTable2);

load(fullfile(pwd, 'FlyPEZ Data\0038000015110252_dataForVisualization.mat')) 
load(fullfile(pwd, 'FlyPEZ Data\0038000015110252_manualAnnotations.mat'))
graphTable1 = join(graphTable,manualAnnotations,'Keys','Row');
load(fullfile(pwd, 'FlyPEZ Data\0038000015110322_dataForVisualization.mat'))
load(fullfile(pwd, 'FlyPEZ Data\0038000015110322_manualAnnotations.mat'))
graphTable2 = join(graphTable,manualAnnotations,'Keys','Row');
GF1DL40 = vertcat(graphTable1,graphTable2);

load(fullfile(pwd, 'FlyPEZ Data\0038000002250253_dataForVisualization.mat'))
load(fullfile(pwd, 'FlyPEZ Data\0038000002250253_manualAnnotations.mat'))
graphTable1 = join(graphTable,manualAnnotations,'Keys','Row');
load(fullfile(pwd, 'FlyPEZ Data\0038000002250321_dataForVisualization.mat'))
load(fullfile(pwd, 'FlyPEZ Data\0038000002250321_manualAnnotations.mat'))
graphTable2 = join(graphTable,manualAnnotations,'Keys','Row');
GF1Gal10 = vertcat(graphTable1,graphTable2);

load(fullfile(pwd, 'FlyPEZ Data\0038000002250252_dataForVisualization.mat'))
load(fullfile(pwd, 'FlyPEZ Data\0038000002250252_manualAnnotations.mat'))
graphTable1 = join(graphTable,manualAnnotations,'Keys','Row');
load(fullfile(pwd, 'FlyPEZ Data\0038000002250322_dataForVisualization.mat'))
load(fullfile(pwd, 'FlyPEZ Data\0038000002250322_manualAnnotations.mat'))
graphTable2 = join(graphTable,manualAnnotations,'Keys','Row');
GF1Gal40 = vertcat(graphTable1,graphTable2);

load(fullfile(pwd, 'FlyPEZ Data\0038000002240253_dataForVisualization.mat'))
load(fullfile(pwd, 'FlyPEZ Data\0038000002240253_manualAnnotations.mat'))
graphTable1 = join(graphTable,manualAnnotations,'Keys','Row');
load(fullfile(pwd, 'FlyPEZ Data\0038000002240321_dataForVisualization.mat'))
load(fullfile(pwd, 'FlyPEZ Data\0038000002240321_manualAnnotations.mat'))
graphTable2 = join(graphTable,manualAnnotations,'Keys','Row');
GF1Kir10 = vertcat(graphTable1,graphTable2);

load(fullfile(pwd, 'FlyPEZ Data\0038000002240252_dataForVisualization.mat'))
load(fullfile(pwd, 'FlyPEZ Data\0038000002240252_manualAnnotations.mat'))
graphTable1 = join(graphTable,manualAnnotations,'Keys','Row');
load(fullfile(pwd, 'FlyPEZ Data\0038000002240322_dataForVisualization.mat'))
load(fullfile(pwd, 'FlyPEZ Data\0038000002240322_manualAnnotations.mat'))
graphTable2 = join(graphTable,manualAnnotations,'Keys','Row');
GF1Kir40 = vertcat(graphTable1,graphTable2);

load(fullfile(pwd, 'FlyPEZ Data\0038000015190253_dataForVisualization.mat'))
load(fullfile(pwd, 'FlyPEZ Data\0038000015190253_manualAnnotations.mat'))
GF2DL10 = join(graphTable,manualAnnotations,'Keys','Row');

load(fullfile(pwd, 'FlyPEZ Data\0038000015190252_dataForVisualization.mat'))
load(fullfile(pwd, 'FlyPEZ Data\0038000015190252_manualAnnotations.mat'))
GF2DL40 = join(graphTable,manualAnnotations,'Keys','Row');

load(fullfile(pwd, 'FlyPEZ Data\0038000015180253_dataForVisualization.mat'))
load(fullfile(pwd, 'FlyPEZ Data\0038000015180253_manualAnnotations.mat'))
GF2Kir10 = join(graphTable,manualAnnotations,'Keys','Row');

load(fullfile(pwd, 'FlyPEZ Data\0038000015180252_dataForVisualization.mat'))
load(fullfile(pwd, 'FlyPEZ Data\0038000015180252_manualAnnotations.mat'))
GF2Kir40 = join(graphTable,manualAnnotations,'Keys','Row');

load(fullfile(pwd, 'FlyPEZ Data\0038000015210253_dataForVisualization.mat'))
load(fullfile(pwd, 'FlyPEZ Data\0038000015210253_manualAnnotations.mat'))
empKir10 = join(graphTable,manualAnnotations,'Keys','Row');

load(fullfile(pwd, 'FlyPEZ Data\0038000015210252_dataForVisualization.mat'))
load(fullfile(pwd, 'FlyPEZ Data\0038000015210252_manualAnnotations.mat'))
empKir40 = join(graphTable,manualAnnotations,'Keys','Row');

GF1X = [1 2 3];
GF2X = [1 3];

GF110 = 100.*[sum(GF1DL10.manualJumpTest)/length(GF1DL10.manualJumpTest) sum(GF1Gal10.manualJumpTest)/length(GF1Gal10.manualJumpTest) sum(GF1Kir10.manualJumpTest)/length(GF1Kir10.manualJumpTest)];
GF140 = 100.*[sum(GF1DL40.manualJumpTest)/length(GF1DL40.manualJumpTest) sum(GF1Gal40.manualJumpTest)/length(GF1Gal40.manualJumpTest) sum(GF1Kir40.manualJumpTest)/length(GF1Kir40.manualJumpTest)];
GF210 = 100.*[sum(GF2DL10.manualJumpTest)/length(GF2DL10.manualJumpTest) sum(GF2Kir10.manualJumpTest)/length(GF2Kir10.manualJumpTest)];
GF240 = 100.*[sum(GF2DL40.manualJumpTest)/length(GF2DL40.manualJumpTest) sum(GF2Kir40.manualJumpTest)/length(GF2Kir40.manualJumpTest)];


figure('Position',[300 300 300 400])
a = plot(GF1X,GF110,'LineWidth',1.5);
a.Color = [0 0 0];
a.Marker = 'o';
a.MarkerFaceColor = [0 0 0];
hold on

c = plot(GF2X,GF210,'LineWidth',1.5);
c.Color = [0 0 0];
c.Marker = 's';
c.MarkerFaceColor = [0 0 0];
hold on

b = plot(GF1X,GF140,'LineWidth',1.5);
b.Color = [0 0 0];
b.Marker = 'o';
b.MarkerFaceColor = [1 1 1];
hold on


d = plot(GF2X,GF240,'LineWidth',1.5);
d.Color = [0 0 0];
d.Marker = 's';
d.MarkerFaceColor = [1 1 1];

ylabel('% Takeoff', 'FontWeight', 'bold')
% xlabel('Peak Acceleration (m/s^{2})', 'FontWeight', 'bold')
legend('GF1>, l/v 10','GF2>, l/v 10','GF1>, l/v 40','GF2>, l/v 40','Location', 'northeast')
% ylim([0 0.6])
% xlim([2.5 4.5])

xticks([1 2 3])
xticklabels(["+","Kir2.1;Gal80ts","Kir2.1"])
xlabel('Effector','FontWeight', 'bold')
ax = gca;
ax.FontSize = 14;

set(gca, 'box', 'off')


% Short Calc
short = 7; %short mode threshold in ms

GF110DL = rmmissing((cell2mat(GF1DL10.frame_of_take_off) - cell2mat(GF1DL10.frame_of_wing_movement)) /6);
GF110DL = 100*(sum(GF110DL < short)/length(GF110DL));

GF110Gal = rmmissing((cell2mat(GF1Gal10.frame_of_take_off) - cell2mat(GF1Gal10.frame_of_wing_movement)) /6);
GF110Gal = 100*(sum(GF110Gal < short)/length(GF110Gal));

GF110Kir = rmmissing((cell2mat(GF1Kir10.frame_of_take_off) - cell2mat(GF1Kir10.frame_of_wing_movement)) /6);
GF110Kir = 100*(sum(GF110Kir < short)/length(GF110Kir));

GF140DL = rmmissing((cell2mat(GF1DL40.frame_of_take_off) - cell2mat(GF1DL40.frame_of_wing_movement)) /6);
GF140DL = 100*(sum(GF140DL < short)/length(GF140DL));

GF140Gal = rmmissing((cell2mat(GF1Gal40.frame_of_take_off) - cell2mat(GF1Gal40.frame_of_wing_movement)) /6);
GF140Gal = 100*(sum(GF140Gal < short)/length(GF140Gal));

GF140Kir = rmmissing((cell2mat(GF1Kir40.frame_of_take_off) - cell2mat(GF1Kir40.frame_of_wing_movement)) /6);
GF140Kir = 100*(sum(GF140Kir < short)/length(GF140Kir));

GF210DL = rmmissing((cell2mat(GF2DL10.frame_of_take_off) - cell2mat(GF2DL10.frame_of_wing_movement)) /6);
GF210DL = 100*(sum(GF210DL < short)/length(GF210DL));

GF210Kir = rmmissing((cell2mat(GF2Kir10.frame_of_take_off) - cell2mat(GF2Kir10.frame_of_wing_movement)) /6);
GF210Kir = 100*(sum(GF210Kir < short)/length(GF210Kir));

GF240DL = rmmissing((cell2mat(GF2DL40.frame_of_take_off) - cell2mat(GF2DL40.frame_of_wing_movement)) /6);
GF240DL = 100*(sum(GF240DL < short)/length(GF240DL));

GF240Kir = rmmissing((cell2mat(GF2Kir40.frame_of_take_off) - cell2mat(GF2Kir40.frame_of_wing_movement)) /6);
GF240Kir = 100*(sum(GF240Kir < short)/length(GF240Kir));

GF110 = [GF110DL GF110Gal GF110Kir];
GF140 = [GF140DL GF140Gal GF140Kir];
GF210 = [GF210DL          GF210Kir];
GF240 = [GF240DL          GF240Kir];

figure('Position',[300 300 300 400])
a = plot(GF1X,GF110,'LineWidth',1.5);
a.Color = [0 0 0];
a.Marker = 'o';
a.MarkerFaceColor = [0 0 0];
hold on

c = plot(GF2X,GF210,'LineWidth',1.5);
c.Color = [0 0 0];
c.Marker = 's';
c.MarkerFaceColor = [0 0 0];
hold on

b = plot(GF1X,GF140,'LineWidth',1.5);
b.Color = [0 0 0];
b.Marker = 'o';
b.MarkerFaceColor = [1 1 1];
hold on


d = plot(GF2X,GF240,'LineWidth',1.5);
d.Color = [0 0 0];
d.Marker = 's';
d.MarkerFaceColor = [1 1 1];

ylabel('% Short Takeoff', 'FontWeight', 'bold')
% xlabel('Peak Acceleration (m/s^{2})', 'FontWeight', 'bold')
%legend('GF1, l/v 10','GF2, l/v 10','GF1, l/v 40','GF2, l/v 40','Location', 'northeast')
% ylim([0 0.6])
% xlim([2.5 4.5])

xticks([1 2 3])
xticklabels(["+","Kir2.1;Gal80ts","Kir2.1"])
xlabel('Effector','FontWeight', 'bold')
ax = gca;
ax.FontSize = 14;

set(gca, 'box', 'off')



%Takeoff time plot

GF110DLtime = rmmissing((cell2mat(GF1DL10.frame_of_take_off) - cell2mat(GF1DL10.frame_of_wing_movement)) /6);
GF140DLtime = rmmissing((cell2mat(GF1DL40.frame_of_take_off) - cell2mat(GF1DL40.frame_of_wing_movement)) /6);
GF210DLtime = rmmissing((cell2mat(GF2DL10.frame_of_take_off) - cell2mat(GF2DL10.frame_of_wing_movement)) /6);
GF240DLtime = rmmissing((cell2mat(GF2DL40.frame_of_take_off) - cell2mat(GF2DL40.frame_of_wing_movement)) /6);
DLplt = [GF110DLtime; GF140DLtime; GF210DLtime ;GF240DLtime];

GF110KIRtime = rmmissing((cell2mat(GF1Kir10.frame_of_take_off) - cell2mat(GF1Kir10.frame_of_wing_movement)) /6);
GF140KIRtime = rmmissing((cell2mat(GF1Kir40.frame_of_take_off) - cell2mat(GF1Kir40.frame_of_wing_movement)) /6);
GF210KIRtime = rmmissing((cell2mat(GF2Kir10.frame_of_take_off) - cell2mat(GF2Kir10.frame_of_wing_movement)) /6);
GF240KIRtime = rmmissing((cell2mat(GF2Kir40.frame_of_take_off) - cell2mat(GF2Kir40.frame_of_wing_movement)) /6);
KIRplt = [GF110KIRtime; GF140KIRtime; GF210KIRtime ;GF240KIRtime];


DL10 = [GF110DLtime; GF210DLtime];
DL40 = [GF140DLtime; GF240DLtime];
Kir10 = [GF110KIRtime; GF210KIRtime];
Kir40 = [GF140KIRtime; GF240KIRtime];

% for i = 1:(length(DL10)-length(Kir10))
%     check = 1;
%     while check == 1
%         ind = round(rand*length(DL10));
%         if ind == 0
%             ind = 1;
%         end
%         if ~isnan(DL10(ind))
%             DL10(ind) = NaN;
%             check = 0;
%         end
%     end
% end
for  i = 1:(length(DL40)-length(DL10))
    check = 1;
    while check == 1
        ind = round(rand*length(DL40));
        if ind == 0
            ind = 1;
        end
        if ~isnan(DL40(ind))
            DL40(ind) = NaN;
            check = 0;
        end
    end
end

for i = 1:(length(Kir40)-length(Kir10))
    check = 1;
    while check == 1
        ind = round(rand*length(Kir40));
        if ind == 0
            ind = 1;
        end
        if ~isnan(Kir40(ind))
            Kir40(ind) = NaN;
            check = 0;
        end
    end
end

% DL10 = rmmissing(DL10);
DL40 = rmmissing(DL40);
Kir40 = rmmissing(Kir40);

DLplt = [DL10;DL40];
KIRplt = [Kir10;Kir40];

edges = 1:2:100;

figure('Position',[300 300 800 250])
% [~,edges] = histcounts(log10(timeplt),50);
% histogram(timeplt,10.^edges)
% set(gca, 'xscale','log')
 histogram(DLplt,edges,'Normalization','probability')
hold on
histogram(KIRplt,edges,'Normalization','probability')
hold on
%plot([7 7],[0 100],'-r','LineWidth',2)
xlim([0,40])
ylim([0,.4])
 yticks([0 .1 .2 .3 .4])

ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
xlabel('Escape Sequence Duration (ms)','FontWeight', 'bold','FontSize',14)
ylabel('Percentage Flies','FontWeight', 'bold','FontSize',14)
legend('GF>+','GF>Kir2.1','Location', 'northeast')
ax = gca;
ax.FontSize = 14;

set(gca, 'box', 'off')


%% Time between stim start and winglift

GF110DLtime = rmmissing((cell2mat(GF1DL10.frame_of_wing_movement) - GF1DL10.stimStart) /6);
GF140DLtime = rmmissing((cell2mat(GF1DL40.frame_of_wing_movement) - GF1DL40.stimStart) /6);
GF210DLtime = rmmissing((cell2mat(GF2DL10.frame_of_wing_movement) - GF2DL10.stimStart) /6);
GF240DLtime = rmmissing((cell2mat(GF2DL40.frame_of_wing_movement) - GF2DL40.stimStart) /6);
DL10plt = [GF110DLtime; GF210DLtime];
DL40plt = [GF140DLtime; GF240DLtime];

DL10plt(DL10plt<=0)=[];
DL40plt(DL40plt<=0)=[];

GF110KIRtime = rmmissing((cell2mat(GF1Kir10.frame_of_wing_movement) - GF1Kir10.stimStart) /6);
GF140KIRtime = rmmissing((cell2mat(GF1Kir40.frame_of_wing_movement) - GF1Kir40.stimStart) /6);
GF210KIRtime = rmmissing((cell2mat(GF2Kir10.frame_of_wing_movement) - GF2Kir10.stimStart) /6);
GF240KIRtime = rmmissing((cell2mat(GF2Kir40.frame_of_wing_movement) - GF2Kir40.stimStart) /6);
KIR10plt = [GF110KIRtime; GF210KIRtime];
KIR40plt = [GF140KIRtime; GF240KIRtime];

KIR10plt(KIR10plt<=0)=[];
KIR40plt(KIR40plt<=0)=[];

emp10plt = rmmissing((cell2mat(empKir10.frame_of_wing_movement) - empKir10.stimStart) /6);
emp40plt = rmmissing((cell2mat(empKir40.frame_of_wing_movement) - empKir40.stimStart) /6);

%emp10plt(emp10plt<=50)=[];
%emp40plt(emp40plt<=100)=[];

edges = 0:4:(ceil(max([DL10plt;KIR10plt]))+4);

figure('Position',[300 300 500 250])
histogram(DL10plt,edges,'Normalization','probability')
hold on
histogram(KIR10plt,edges,'Normalization','probability')
hold on
%plot([mean(DL10plt) mean(DL10plt)],[0,1],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
%plot([mean(DL10plt)-std(DL10plt) mean(DL10plt)+std(DL10plt)],[.2,.2],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
%plot([mean(KIR10plt) mean(KIR10plt)],[0,1],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
hold on 
%plot([mean(KIR10plt)-std(KIR10plt) mean(KIR10plt)+std(KIR10plt)],[.21,.21],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
xlim([0,150])
xlim([0,150])
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
ylim([0,.28])

xlabel('Takeoff Latency Using Start of Wing Lift (ms)','FontWeight', 'bold','FontSize',14)
ylabel('Percentage Flies','FontWeight', 'bold','FontSize',14)
legend('GF>+, n=80','GF>Kir2.1, n=29','Location', 'northwest')
title('l/v = 10')
ax = gca;
ax.FontSize = 14;

set(gca, 'box', 'off')

% Box plot for lv 10 winglift
g1 = repmat({'DL'},length(DL10plt),1);
g2 = repmat({'KIR'},length(KIR10plt),1);
g = [g2; g1];
figure('Position',[300 300 500 150])
boxplot([KIR10plt;DL10plt],g,'Orientation','horizontal')
xlim([0 150])
set(gca, 'box', 'off')

edges = 0:13:(ceil(max([DL40plt;KIR40plt]))+13);
figure('Position',[300 300 500 250])
histogram(DL40plt,edges,'Normalization','probability')
hold on
histogram(KIR40plt,edges,'Normalization','probability')
hold on
%plot([mean(DL40plt) mean(DL40plt)],[0,1],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
%plot([mean(DL40plt)-std(DL40plt) mean(DL40plt)+std(DL40plt)],[.12,.12],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
%plot([mean(KIR40plt) mean(KIR40plt)],[0,1],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
hold on 
%plot([mean(KIR40plt)-std(KIR40plt) mean(KIR40plt)+std(KIR40plt)],[.13,.13],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
xlim([0,500])
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
ylim([0,.2])

xlabel('Takeoff Latency Using Start of Wing Lift (ms)','FontWeight', 'bold','FontSize',14)
ylabel('Percentage Flies','FontWeight', 'bold','FontSize',14)
legend('GF>+, n=141','GF>Kir2.1, n=48','Location', 'northwest')
title('l/v = 40')
ax = gca;
ax.FontSize = 14;

set(gca, 'box', 'off')

% Box plot for lv 40 winglift
g1 = repmat({'DL'},length(DL40plt),1);
g2 = repmat({'KIR'},length(KIR40plt),1);
g = [g2; g1];
figure('Position',[300 300 500 150])
boxplot([KIR40plt;DL40plt],g,'Orientation','horizontal')
xlim([0 500])
set(gca, 'box', 'off')

% Empty x Kir vs GF x Kir

edges = 0:4:(ceil(max([emp10plt;KIR10plt]))+4);

figure('Position',[300 300 500 250])
histogram(emp10plt,edges,'Normalization','probability')
hold on
histogram(KIR10plt,edges,'Normalization','probability')
hold on
plot([mean(emp10plt) mean(emp10plt)],[0,1],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
plot([mean(emp10plt)-std(emp10plt) mean(emp10plt)+std(emp10plt)],[.2,.2],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
plot([mean(KIR10plt) mean(KIR10plt)],[0,1],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
hold on 
plot([mean(KIR10plt)-std(KIR10plt) mean(KIR10plt)+std(KIR10plt)],[.21,.21],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
xlim([0,150])
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
ylim([0,.22])

xlabel('Takeoff Latency Using Start of Wing Lift (ms)','FontWeight', 'bold','FontSize',14)
ylabel('Percentage Flies','FontWeight', 'bold','FontSize',14)
legend('Empty>Kir2.1, n=20','GF>Kir2.1, n=29','Location', 'northwest')
title('l/v = 10')
ax = gca;
ax.FontSize = 14;

set(gca, 'box', 'off')

edges = 0:13:(ceil(max([emp40plt;KIR40plt]))+13);
figure('Position',[300 300 500 250])
histogram(emp40plt,edges,'Normalization','probability')
hold on
histogram(KIR40plt,edges,'Normalization','probability')
hold on
plot([mean(emp40plt) mean(emp40plt)],[0,1],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
plot([mean(emp40plt)-std(emp40plt) mean(emp40plt)+std(emp40plt)],[.12,.12],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
plot([mean(KIR40plt) mean(KIR40plt)],[0,1],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
hold on 
plot([mean(KIR40plt)-std(KIR40plt) mean(KIR40plt)+std(KIR40plt)],[.13,.13],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
xlim([0,500])
yticks([0 .05 .1 .15])
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
ylim([0,.14])

xlabel('Takeoff Latency Using Start of Wing Lift (ms)','FontWeight', 'bold','FontSize',14)
ylabel('Percentage Flies','FontWeight', 'bold','FontSize',14)
legend('Empty>Kir2.1, n=40','GF>Kir2.1, n=48','Location', 'northwest')
title('l/v = 40')
ax = gca;
ax.FontSize = 14;

set(gca, 'box', 'off')


%% Time between stim start and takeoff

GF110DLtime = rmmissing((cell2mat(GF1DL10.frame_of_take_off) - GF1DL10.stimStart) /6);
GF140DLtime = rmmissing((cell2mat(GF1DL40.frame_of_take_off) - GF1DL40.stimStart) /6);
GF210DLtime = rmmissing((cell2mat(GF2DL10.frame_of_take_off) - GF2DL10.stimStart) /6);
GF240DLtime = rmmissing((cell2mat(GF2DL40.frame_of_take_off) - GF2DL40.stimStart) /6);
DL10plt = [GF110DLtime; GF210DLtime];
DL40plt = [GF140DLtime; GF240DLtime];

DL10plt(DL10plt<=0)=[];
DL40plt(DL40plt<=0)=[];

GF110KIRtime = rmmissing((cell2mat(GF1Kir10.frame_of_take_off) - GF1Kir10.stimStart) /6);
GF140KIRtime = rmmissing((cell2mat(GF1Kir40.frame_of_take_off) - GF1Kir40.stimStart) /6);
GF210KIRtime = rmmissing((cell2mat(GF2Kir10.frame_of_take_off) - GF2Kir10.stimStart) /6);
GF240KIRtime = rmmissing((cell2mat(GF2Kir40.frame_of_take_off) - GF2Kir40.stimStart) /6);
KIR10plt = [GF110KIRtime; GF210KIRtime];
KIR40plt = [GF140KIRtime; GF240KIRtime];

KIR10plt(KIR10plt<=0)=[];
KIR40plt(KIR40plt<=0)=[];

emp10plt = rmmissing((cell2mat(empKir10.frame_of_take_off) - empKir10.stimStart) /6);
emp40plt = rmmissing((cell2mat(empKir40.frame_of_take_off) - empKir40.stimStart) /6);

%emp10plt(emp10plt<=50)=[];
%emp40plt(emp40plt<=100) = [];

edges = 0:4:(ceil(max([DL10plt;KIR10plt]))+4);

figure('Position',[300 300 500 250])
histogram(DL10plt,edges,'Normalization','probability')
hold on
histogram(KIR10plt,edges,'Normalization','probability')
hold on
%plot([mean(DL10plt) mean(DL10plt)],[0,1],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
%plot([mean(DL10plt)-std(DL10plt) mean(DL10plt)+std(DL10plt)],[.2,.2],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
%plot([mean(KIR10plt) mean(KIR10plt)],[0,1],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
hold on 
%plot([mean(KIR10plt)-std(KIR10plt) mean(KIR10plt)+std(KIR10plt)],[.21,.21],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
xlim([0,150])
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
ylim([0,.22])

xlabel('Takeoff Latency (ms)','FontWeight', 'bold','FontSize',14)
ylabel('Percentage Flies','FontWeight', 'bold','FontSize',14)
legend('GF>DL, n=80','GF>Kir2.1, n=29','Location', 'northwest')
title('l/v = 10')
ax = gca;
ax.FontSize = 14;

set(gca, 'box', 'off')

% Box plot for lv 10 TO
g1 = repmat({'DL'},length(DL10plt),1);
g2 = repmat({'KIR'},length(KIR10plt),1);
g = [g2; g1];
figure('Position',[300 300 500 150])
boxplot([KIR10plt;DL10plt],g,'Orientation','horizontal')
xlim([0 150])
set(gca, 'box', 'off')

edges = 0:13:(ceil(max([DL40plt;KIR40plt]))+13);
figure('Position',[300 300 500 250])
histogram(DL40plt,edges,'Normalization','probability')
hold on
histogram(KIR40plt,edges,'Normalization','probability')
hold on
%plot([mean(DL40plt) mean(DL40plt)],[0,1],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
%plot([mean(DL40plt)-std(DL40plt) mean(DL40plt)+std(DL40plt)],[.12,.12],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
%plot([mean(KIR40plt) mean(KIR40plt)],[0,1],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
hold on 
%plot([mean(KIR40plt)-std(KIR40plt) mean(KIR40plt)+std(KIR40plt)],[.13,.13],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
xlim([0,500])
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
ylim([0,.14])

xlabel('Takeoff Latency (ms)','FontWeight', 'bold','FontSize',14)
ylabel('Percentage Flies','FontWeight', 'bold','FontSize',14)
legend('GF>DL, n=141','GF>Kir2.1, n=48','Location', 'northwest')
title('l/v = 40')
ax = gca;
ax.FontSize = 14;

set(gca, 'box', 'off')

% Box plot for lv 40 TO
g1 = repmat({'DL'},length(DL40plt),1);
g2 = repmat({'KIR'},length(KIR40plt),1);
g = [g2; g1];
figure('Position',[300 300 500 150])
boxplot([KIR40plt;DL40plt],g,'Orientation','horizontal')
xlim([0 500])
set(gca, 'box', 'off')

figure('Position',[300 300 500 150])
plot([mean(KIR40plt)-std(KIR40plt) mean(KIR40plt)+std(KIR40plt)],[1,1],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
hold on
plot([mean(KIR40plt) mean(KIR40plt)],[0.75,1.25],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
hold on
plot([mean(DL40plt)-std(DL40plt) mean(DL40plt)+std(DL40plt)],[2,2],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
plot([mean(DL40plt) mean(DL40plt)],[1.75,2.25],'-','LineWidth',2,'Color',[0 0.4470 0.7410])

xlabel('time (ms)')
set(gca, 'box', 'off')
ylim([0,3])
xlim([0,500])

figure('Position',[300 300 500 150])
plot([mean(KIR10plt)-std(KIR10plt) mean(KIR10plt)+std(KIR10plt)],[1,1],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
hold on
plot([mean(KIR10plt) mean(KIR10plt)],[0.75,1.25],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
hold on
plot([mean(DL10plt)-std(DL10plt) mean(DL10plt)+std(DL10plt)],[2,2],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
plot([mean(DL10plt) mean(DL10plt)],[1.75,2.25],'-','LineWidth',2,'Color',[0 0.4470 0.7410])

xlabel('time (ms)')
set(gca, 'box', 'off')
ylim([0,3])
xlim([0,150])

% Empty x Kir vs GF x Kir

edges = 0:4:ceil(max([emp10plt;KIR10plt]));

figure('Position',[300 300 500 250])
histogram(emp10plt,edges,'Normalization','probability')
hold on
histogram(KIR10plt,edges,'Normalization','probability')
hold on
plot([mean(emp10plt) mean(emp10plt)],[0,1],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
plot([mean(emp10plt)-std(emp10plt) mean(emp10plt)+std(emp10plt)],[.2,.2],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
plot([mean(KIR10plt) mean(KIR10plt)],[0,1],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
hold on 
plot([mean(KIR10plt)-std(KIR10plt) mean(KIR10plt)+std(KIR10plt)],[.21,.21],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
xlim([0,150])
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
ylim([0,.22])

xlabel('Takeoff Latency (ms)','FontWeight', 'bold','FontSize',14)
ylabel('Percentage Flies','FontWeight', 'bold','FontSize',14)
legend('Empty>Kir2.1, n=20','GF>Kir2.1, n=29','Location', 'northwest')
title('l/v = 10')
ax = gca;
ax.FontSize = 14;

set(gca, 'box', 'off')

edges = 0:13:ceil(max([emp40plt;KIR40plt]));
figure('Position',[300 300 500 250])
histogram(emp40plt,edges,'Normalization','probability')
hold on
histogram(KIR40plt,edges,'Normalization','probability')
hold on
plot([mean(emp40plt) mean(emp40plt)],[0,1],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
plot([mean(emp40plt)-std(emp40plt) mean(emp40plt)+std(emp40plt)],[.12,.12],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
plot([mean(KIR40plt) mean(KIR40plt)],[0,1],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
hold on 
plot([mean(KIR40plt)-std(KIR40plt) mean(KIR40plt)+std(KIR40plt)],[.13,.13],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
xlim([0,500])
yticks([0 .05 .1 .15])
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
ylim([0,.14])

xlabel('Takeoff Latency(ms)','FontWeight', 'bold','FontSize',14)
ylabel('Percentage Flies','FontWeight', 'bold','FontSize',14)
legend('Empty>Kir2.1, n=40','GF>Kir2.1, n=48','Location', 'northwest')
title('l/v = 40')
ax = gca;
ax.FontSize = 14;

set(gca, 'box', 'off')


%% Time between stim start and legpush

GF110DLtime = rmmissing((cell2mat(GF1DL10.frame_of_leg_push) - GF1DL10.stimStart) /6);
GF140DLtime = rmmissing((cell2mat(GF1DL40.frame_of_leg_push) - GF1DL40.stimStart) /6);
GF210DLtime = rmmissing((cell2mat(GF2DL10.frame_of_leg_push) - GF2DL10.stimStart) /6);
GF240DLtime = rmmissing((cell2mat(GF2DL40.frame_of_leg_push) - GF2DL40.stimStart) /6);
DL10plt = [GF110DLtime; GF210DLtime];
DL40plt = [GF140DLtime; GF240DLtime];

DL10plt(DL10plt<=0)=[];
DL40plt(DL40plt<=0)=[];

GF110KIRtime = rmmissing((cell2mat(GF1Kir10.frame_of_leg_push) - GF1Kir10.stimStart) /6);
GF140KIRtime = rmmissing((cell2mat(GF1Kir40.frame_of_leg_push) - GF1Kir40.stimStart) /6);
GF210KIRtime = rmmissing((cell2mat(GF2Kir10.frame_of_leg_push) - GF2Kir10.stimStart) /6);
GF240KIRtime = rmmissing((cell2mat(GF2Kir40.frame_of_leg_push) - GF2Kir40.stimStart) /6);
KIR10plt = [GF110KIRtime; GF210KIRtime];
KIR40plt = [GF140KIRtime; GF240KIRtime];

KIR10plt(KIR10plt<=0)=[];
KIR40plt(KIR40plt<=0)=[];

emp10plt = rmmissing((cell2mat(empKir10.frame_of_leg_push) - empKir10.stimStart) /6);
emp40plt = rmmissing((cell2mat(empKir40.frame_of_leg_push) - empKir40.stimStart) /6);

%emp10plt(emp10plt<=50)=[];
%emp40plt(emp40plt<=100)=[];

edges = 0:4:(ceil(max([DL10plt;KIR10plt]))+4);

figure('Position',[300 300 500 250])
histogram(DL10plt,edges,'Normalization','probability')
hold on
histogram(KIR10plt,edges,'Normalization','probability')
hold on
%plot([mean(DL10plt) mean(DL10plt)],[0,1],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
%plot([mean(DL10plt)-std(DL10plt) mean(DL10plt)+std(DL10plt)],[.2,.2],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
%plot([mean(KIR10plt) mean(KIR10plt)],[0,1],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
hold on 
%plot([mean(KIR10plt)-std(KIR10plt) mean(KIR10plt)+std(KIR10plt)],[.21,.21],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
xlim([0,150])
xlim([0,150])
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
ylim([0,.28])

xlabel('Leg Push Latency (ms)','FontWeight', 'bold','FontSize',14)
ylabel('Percentage Flies','FontWeight', 'bold','FontSize',14)
legend('GF>+, n=80','GF>Kir2.1, n=29','Location', 'northwest')
title('l/v = 10')
ax = gca;
ax.FontSize = 14;

set(gca, 'box', 'off')

% Box plot for lv 10 winglift
g1 = repmat({'DL'},length(DL10plt),1);
g2 = repmat({'KIR'},length(KIR10plt),1);
g = [g2; g1];
figure('Position',[300 300 500 150])
boxplot([KIR10plt;DL10plt],g,'Orientation','horizontal')
xlim([0 150])
set(gca, 'box', 'off')

edges = 0:13:(ceil(max([DL40plt;KIR40plt]))+13);
figure('Position',[300 300 500 250])
histogram(DL40plt,edges,'Normalization','probability')
hold on
histogram(KIR40plt,edges,'Normalization','probability')
hold on
%plot([mean(DL40plt) mean(DL40plt)],[0,1],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
%plot([mean(DL40plt)-std(DL40plt) mean(DL40plt)+std(DL40plt)],[.12,.12],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
%plot([mean(KIR40plt) mean(KIR40plt)],[0,1],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
hold on 
%plot([mean(KIR40plt)-std(KIR40plt) mean(KIR40plt)+std(KIR40plt)],[.13,.13],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
xlim([0,500])
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
ylim([0,.2])

xlabel('Leg Push Latency (ms)','FontWeight', 'bold','FontSize',14)
ylabel('Percentage Flies','FontWeight', 'bold','FontSize',14)
legend('GF>+, n=141','GF>Kir2.1, n=48','Location', 'northwest')
title('l/v = 40')
ax = gca;
ax.FontSize = 14;

set(gca, 'box', 'off')

% Box plot for lv 40 winglift
g1 = repmat({'DL'},length(DL40plt),1);
g2 = repmat({'KIR'},length(KIR40plt),1);
g = [g2; g1];
figure('Position',[300 300 500 150])
boxplot([KIR40plt;DL40plt],g,'Orientation','horizontal')
xlim([0 500])
set(gca, 'box', 'off')

% Empty x Kir vs GF x Kir

edges = 0:4:(ceil(max([emp10plt;KIR10plt]))+4);

figure('Position',[300 300 500 250])
histogram(emp10plt,edges,'Normalization','probability')
hold on
histogram(KIR10plt,edges,'Normalization','probability')
hold on
plot([mean(emp10plt) mean(emp10plt)],[0,1],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
plot([mean(emp10plt)-std(emp10plt) mean(emp10plt)+std(emp10plt)],[.2,.2],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
plot([mean(KIR10plt) mean(KIR10plt)],[0,1],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
hold on 
plot([mean(KIR10plt)-std(KIR10plt) mean(KIR10plt)+std(KIR10plt)],[.21,.21],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
xlim([0,150])
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
ylim([0,.22])

xlabel('Leg Push Latency (ms)','FontWeight', 'bold','FontSize',14)
ylabel('Percentage Flies','FontWeight', 'bold','FontSize',14)
legend('Empty>Kir2.1, n=20','GF>Kir2.1, n=29','Location', 'northwest')
title('l/v = 10')
ax = gca;
ax.FontSize = 14;

set(gca, 'box', 'off')

edges = 0:13:(ceil(max([emp40plt;KIR40plt]))+13);
figure('Position',[300 300 500 250])
histogram(emp40plt,edges,'Normalization','probability')
hold on
histogram(KIR40plt,edges,'Normalization','probability')
hold on
plot([mean(emp40plt) mean(emp40plt)],[0,1],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
plot([mean(emp40plt)-std(emp40plt) mean(emp40plt)+std(emp40plt)],[.12,.12],'-','LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
plot([mean(KIR40plt) mean(KIR40plt)],[0,1],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
hold on 
plot([mean(KIR40plt)-std(KIR40plt) mean(KIR40plt)+std(KIR40plt)],[.13,.13],'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
xlim([0,500])
yticks([0 .05 .1 .15])
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
ylim([0,.14])

xlabel('Leg Push Latency(ms)','FontWeight', 'bold','FontSize',14)
ylabel('Percentage Flies','FontWeight', 'bold','FontSize',14)
legend('Empty>Kir2.1, n=40','GF>Kir2.1, n=48','Location', 'northwest')
title('l/v = 40')
ax = gca;
ax.FontSize = 14;

set(gca, 'box', 'off')



%% STEP 29:Individual distance plots

%Setup speed variables for plotting
distp1 = padcat(flipud(Type(1).Data(1).distVec_Dams), flipud(Type(1).Data(2).distVec_Dams), flipud(Type(1).Data(3).distVec_Dams), flipud(Type(1).Data(4).distVec_Dams), flipud(Type(1).Data(5).distVec_Dams)...
    ,flipud(Type(1).Data(6).distVec_Dams), flipud(Type(1).Data(7).distVec_Dams), flipud(Type(1).Data(8).distVec_Dams), flipud(Type(1).Data(6).distVec_Dams), flipud(Type(1).Data(10).distVec_Dams), flipud(Type(1).Data(11).distVec_Dams)...
    ,flipud(Type(1).Data(12).distVec_Dams), flipud(Type(1).Data(13).distVec_Dams));
distp1 = flipud(distp1); %concat dist matricies
numnansp1 = sum(~isnan(distp1),2); %number of traces averaged as a function of time
meanp1Vec = nanmean(distp1, 2); %mean dist as a function of time
stdp1Vec = nanstd(distp1, 0, 2); %standard deviation of dist

distp2 = padcat(flipud(Type(2).Data(1).distVec_Dams), flipud(Type(2).Data(2).distVec_Dams), flipud(Type(2).Data(3).distVec_Dams), flipud(Type(2).Data(4).distVec_Dams), flipud(Type(2).Data(5).distVec_Dams)...
    , flipud(Type(2).Data(6).distVec_Dams));
distp2 = flipud(distp2);
numnansp2 = sum(~isnan(distp2),2);
meanp2Vec = nanmean(distp2, 2);
stdp2Vec = nanstd(distp2, 0, 2);

distp3 = padcat(flipud(Type(3).Data(1).distVec_Dams), flipud(Type(3).Data(2).distVec_Dams), flipud(Type(3).Data(3).distVec_Dams), flipud(Type(3).Data(4).distVec_Dams), flipud(Type(3).Data(5).distVec_Dams));
distp3 = flipud(distp3);
numnansp3 = sum(~isnan(distp3),2);
meanp3Vec = nanmean(distp3, 2);
stdp3Vec = nanstd(distp3, 0, 2);

distp4 = padcat(flipud(Type(4).Data(1).distVec_Dams), flipud(Type(4).Data(2).distVec_Dams), flipud(Type(4).Data(3).distVec_Dams), flipud(Type(4).Data(4).distVec_Dams), flipud(Type(4).Data(5).distVec_Dams)...
    ,flipud(Type(4).Data(6).distVec_Dams), flipud(Type(4).Data(7).distVec_Dams));
distp4 = flipud(distp4);
numnansp4 = sum(~isnan(distp4),2);
meanp4Vec = nanmean(distp4, 2);
stdp4Vec = nanstd(distp4, 0, 2);

linecolor = [[1 0 0] ;[0 0 1] ;[1 .5 0] ;[0 1 0]];

% Plot DL - Capture dist
figure('Position',[500 300 500 500])

% Uncomment for raw traces
for d = 1:length(Type(3).Data)
    
    plot((Type(3).Data(d).t - Type(3).Data(d).t(end)) * 1000,Type(3).Data(d).distVec_Dams(1:length(Type(3).Data(d).t))*1000, 'Color',linecolor(3,:));
    hold on

timedata = (Type(3).Data(d).t - Type(3).Data(d).t(end)) * 1000;
distdata = rad2deg(Type(3).Data(d).damAzi(1:length(Type(3).Data(d).t)));
combined = [timedata distdata];
a = 0;
end

shadedErrorBar(((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000), meanp3Vec(1:length(Type(3).Data(1).t))*1000, stdp3Vec(1:length(Type(3).Data(1).t))*1000,'lineprops',{'Color',linecolor(3,:),'LineWidth',2},'transparent', 1); %shaded area for standard dev


%title('GF1>+ - Capture','FontSize',14)
ylabel('Distance from Fly (mm)', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[0,50],'Xlim',[-100 0])
set(gca, 'box', 'off')
%ax1 = gca;                   % gca = get current axis
%ax1.YAxis.Visible = 'off';   % remove y-axis
%ax1.XAxis.Visible = 'off';
hold off

% Plot Kir - Capture dist
figure('Position',[500 300 500 500])

% Uncomment for raw traces
for d = 1:length(Type(1).Data)
    
    plot((Type(1).Data(d).t - Type(1).Data(d).t(end)) * 1000,Type(1).Data(d).distVec_Dams(1:length(Type(1).Data(d).t))*1000, 'Color',linecolor(1,:));
    hold on


timedata = (Type(1).Data(d).t - Type(1).Data(d).t(end)) * 1000;
distdata = rad2deg(Type(1).Data(d).damAzi(1:length(Type(1).Data(d).t)));
combined = [timedata distdata];
a = 0;
end

shadedErrorBar(((Type(1).Data(3).t - Type(1).Data(3).t(end)) * 1000), meanp1Vec(1:length(Type(1).Data(3).t))*1000, stdp1Vec(1:length(Type(1).Data(3).t))*1000,'lineprops',{'Color',linecolor(1,:),'LineWidth',2},'transparent', 1); %shaded area for standard dev


%title('GF1>Kir2.1 - Capture','FontSize',14)
ylabel('Distance from Fly (mm)', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[0,50],'Xlim',[-100 0])
set(gca, 'box', 'off')

ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';

% Plot DL - Escape dist
figure('Position',[500 300 500 500])

% Uncomment for raw traces
for d = 1:length(Type(4).Data)
    
    plot((Type(4).Data(d).t - Type(4).Data(d).t(end)) * 1000,Type(4).Data(d).distVec_Dams(1:length(Type(4).Data(d).t))*1000, 'Color',linecolor(4,:));
    hold on

    timedata = (Type(4).Data(d).t - Type(4).Data(d).t(end)) * 1000;
    distdata = rad2deg(Type(4).Data(d).damAzi(1:length(Type(4).Data(d).t)));
    combined = [timedata distdata];
    a = 0;
end

shadedErrorBar(((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000), meanp4Vec(1:length(Type(4).Data(6).t))*1000, stdp4Vec(1:length(Type(4).Data(6).t))*1000,'lineprops',{'Color',linecolor(4,:),'LineWidth',2},'transparent', 1); %shaded area for standard dev


%title('GF1>+ - Escape','FontSize',14)
ylabel('Distance from Fly (mm)', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[0,50],'Xlim',[-100 0])
set(gca, 'box', 'off')
ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';
hold off

% Plot Kir - Escape dist
figure('Position',[500 300 500 500])

% Uncomment for raw traces
for d = 1:length(Type(2).Data)
    
    plot((Type(2).Data(d).t - Type(2).Data(d).t(end)) * 1000,Type(2).Data(d).distVec_Dams(1:length(Type(2).Data(d).t))*1000, 'Color',linecolor(2,:));
    hold on

      timedata = (Type(2).Data(d).t - Type(2).Data(d).t(end)) * 1000;
    distdata = rad2deg(Type(2).Data(d).damAzi(1:length(Type(2).Data(d).t)));
    combined = [timedata distdata];
    a = 0;
end

shadedErrorBar(((Type(2).Data(5).t - Type(2).Data(5).t(end)) * 1000), meanp2Vec(1:length(Type(2).Data(5).t))*1000, stdp2Vec(1:length(Type(2).Data(5).t))*1000,'lineprops',{'Color',linecolor(2,:),'LineWidth',2},'transparent', 1); %shaded area for standard dev


%title('GF1>Kir2.1 - Escape','FontSize',14)
ylabel('Distance from Fly (mm)', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[0,50],'Xlim',[-100 0])
set(gca, 'box', 'off')
ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';


%% STEP 30:Individual speed plots

%Setup speed variables for plotting
speedp1 = padcat(flipud(Type(1).Data(1).speedVec_Dams), flipud(Type(1).Data(2).speedVec_Dams), flipud(Type(1).Data(3).speedVec_Dams), flipud(Type(1).Data(4).speedVec_Dams), flipud(Type(1).Data(5).speedVec_Dams)...
    ,flipud(Type(1).Data(6).speedVec_Dams), flipud(Type(1).Data(7).speedVec_Dams), flipud(Type(1).Data(8).speedVec_Dams), flipud(Type(1).Data(6).speedVec_Dams), flipud(Type(1).Data(10).speedVec_Dams), flipud(Type(1).Data(11).speedVec_Dams)...
    ,flipud(Type(1).Data(12).speedVec_Dams), flipud(Type(1).Data(13).speedVec_Dams));
speedp1 = flipud(speedp1); %concat speed matricies
numnansp1 = sum(~isnan(speedp1),2); %number of traces averaged as a function of time
meanp1Vec = nanmean(speedp1, 2); %mean speed as a function of time
stdp1Vec = nanstd(speedp1, 0, 2); %standard deviation of speed

indlow = find(maxSpeed{1}<mean(maxSpeed{1}));
indhig = find(maxSpeed{1}>mean(maxSpeed{1}));

speedl = speedp1(:,indlow);
meanlVec = nanmean(speedl,2);
stdlVec = nanstd(speedl,0,2);

speedh = speedp1(:,indhig);
meanhVec = nanmean(speedh,2);
stdhVec = nanstd(speedh,0,2);

speedp2 = padcat(flipud(Type(2).Data(1).speedVec_Dams), flipud(Type(2).Data(2).speedVec_Dams), flipud(Type(2).Data(3).speedVec_Dams), flipud(Type(2).Data(4).speedVec_Dams), flipud(Type(2).Data(5).speedVec_Dams)...
    , flipud(Type(2).Data(6).speedVec_Dams));
speedp2 = flipud(speedp2);
numnansp2 = sum(~isnan(speedp2),2);
meanp2Vec = nanmean(speedp2, 2);
stdp2Vec = nanstd(speedp2, 0, 2);

speedp3 = padcat(flipud(Type(3).Data(1).speedVec_Dams), flipud(Type(3).Data(2).speedVec_Dams), flipud(Type(3).Data(3).speedVec_Dams), flipud(Type(3).Data(4).speedVec_Dams), flipud(Type(3).Data(5).speedVec_Dams));
speedp3 = flipud(speedp3);
numnansp3 = sum(~isnan(speedp3),2);
meanp3Vec = nanmean(speedp3, 2);
stdp3Vec = nanstd(speedp3, 0, 2);

speedp4 = padcat(flipud(Type(4).Data(1).speedVec_Dams), flipud(Type(4).Data(2).speedVec_Dams), flipud(Type(4).Data(3).speedVec_Dams), flipud(Type(4).Data(4).speedVec_Dams), flipud(Type(4).Data(5).speedVec_Dams)...
    ,flipud(Type(4).Data(6).speedVec_Dams), flipud(Type(4).Data(7).speedVec_Dams));
speedp4 = flipud(speedp4);
numnansp4 = sum(~isnan(speedp4),2);
meanp4Vec = nanmean(speedp4, 2);
stdp4Vec = nanstd(speedp4, 0, 2);

linecolor = [[1 0 0] ;[0 0 1] ;[1 .5 0] ;[0 1 0]];

% Plot DL - Capture Speed
figure('Position',[500 300 500 500])

% Uncomment for raw traces
for d = 1:length(Type(3).Data)
    
    plot((Type(3).Data(d).t - Type(3).Data(d).t(end)) * 1000,Type(3).Data(d).speedVec_Dams(1:length(Type(3).Data(d).t)), 'Color',linecolor(3,:));
    hold on
end

shadedErrorBar(((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000), meanp3Vec(1:length(Type(3).Data(1).t)), stdp3Vec(1:length(Type(3).Data(1).t)),'lineprops',{'Color',linecolor(3,:),'LineWidth',2},'transparent', 1); %shaded area for standard dev

for d = 1:length(Type(3).Data)
    time = (Type(3).Data(d).t - Type(3).Data(d).t(end)) * 1000;
    speed = Type(3).Data(d).speedVec_Dams(1:length(Type(3).Data(d).t));
    [m,ind] = max(speed);
    maxtime(d) = time(ind);
    plot((Type(3).Data(d).t - Type(3).Data(d).t(end)) * 1000,Type(3).Data(d).speedVec_Dams(1:length(Type(3).Data(d).t)), 'Color',linecolor(3,:));
    hold on
end


%title('GF1>+ - Capture','FontSize',14)
ylabel('Speed (m/s)', 'FontWeight', 'bold','FontSize',14)
xlabel('time(ms)', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[0,0.325],'Xlim',[-100 0])
hold on 
plot([-29,-29],[0,.3],'-k')
yticks([0 .1 .2 .3 .4])
set(gca, 'box', 'off')
%ax1 = gca;                   % gca = get current axis
%ax1.YAxis.Visible = 'off';   % remove y-axis
%ax1.XAxis.Visible = 'off';
hold off

% Plot Kir - Capture Speed
figure('Position',[500 300 500 500])

% Uncomment for raw traces
for d = 1:length(Type(1).Data)
    
    plot((Type(1).Data(d).t - Type(1).Data(d).t(end)) * 1000,Type(1).Data(d).speedVec_Dams(1:length(Type(1).Data(d).t)), 'Color',linecolor(1,:));
    hold on
end

shadedErrorBar(((Type(1).Data(3).t - Type(1).Data(3).t(end)) * 1000), meanp1Vec(1:length(Type(1).Data(3).t)), stdp1Vec(1:length(Type(1).Data(3).t)),'lineprops',{'Color',linecolor(1,:),'LineWidth',2},'transparent', 1); %shaded area for standard dev


%title('GF1>Kir2.1 - Capture','FontSize',14)
ylabel('Speed (m/s)', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[0,0.325],'Xlim',[-100 0])
set(gca, 'box', 'off')
ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';
yticks([0 .1 .2 .3 .4])

% Plot DL - Escape Speed
figure('Position',[500 300 500 500])

% Uncomment for raw traces
for d = 1:length(Type(4).Data)
    
    plot((Type(4).Data(d).t - Type(4).Data(d).t(end)) * 1000,Type(4).Data(d).speedVec_Dams(1:length(Type(4).Data(d).t)), 'Color',linecolor(4,:));
    hold on
end

shadedErrorBar(((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000), meanp4Vec(1:length(Type(4).Data(6).t)), stdp4Vec(1:length(Type(4).Data(6).t)),'lineprops',{'Color',linecolor(4,:),'LineWidth',2},'transparent', 1); %shaded area for standard dev


%title('GF1>+ - Escape','FontSize',14)
ylabel('Speed (m/s)', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[0,0.325],'Xlim',[-100 0])
set(gca, 'box', 'off')
ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';
yticks([0 .1 .2 .3 .4])
hold off

% Plot Kir - Escape Speed
figure('Position',[500 300 500 500])

% Uncomment for raw traces
for d = 1:length(Type(2).Data)
    
    plot((Type(2).Data(d).t - Type(2).Data(d).t(end)) * 1000,Type(2).Data(d).speedVec_Dams(1:length(Type(2).Data(d).t)), 'Color',linecolor(2,:));
    hold on
end

shadedErrorBar(((Type(2).Data(5).t - Type(2).Data(5).t(end)) * 1000), meanp2Vec(1:length(Type(2).Data(5).t)), stdp2Vec(1:length(Type(2).Data(5).t)),'lineprops',{'Color',linecolor(2,:),'LineWidth',2},'transparent', 1); %shaded area for standard dev


%title('GF1>Kir2.1 - Escape','FontSize',14)
ylabel('Speed (m/s)', 'FontWeight', 'bold','FontSize',14)
xlabel('Time (ms)', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[0,0.325],'Xlim',[-100 0])
set(gca, 'box', 'off')
ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';
yticks([0 .1 .2 .3 .4])


%Split kir cap into high and low velocity traces

% Plot Kir - Capture Speed High
figure('Position',[500 300 500 500])

% Uncomment for raw traces
for d = 1:length(indhig)
    
    plot((Type(1).Data(indhig(d)).t - Type(1).Data(indhig(d)).t(end)) * 1000,Type(1).Data(indhig(d)).speedVec_Dams(1:length(Type(1).Data(indhig(d)).t)), 'Color',linecolor(1,:));
    hold on
end

shadedErrorBar(((Type(1).Data(3).t - Type(1).Data(3).t(end)) * 1000), meanhVec(1:length(Type(1).Data(3).t)), stdhVec(1:length(Type(1).Data(3).t)),'lineprops',{'Color',linecolor(1,:),'LineWidth',2},'transparent', 1); %shaded area for standard dev


%title('GF1>Kir2.1 - Capture','FontSize',14)
ylabel('Speed (m/s)', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[0,0.325],'Xlim',[-100 0])
set(gca, 'box', 'off')
%ax1 = gca;                   % gca = get current axis
%ax1.YAxis.Visible = 'off';   % remove y-axis
%ax1.XAxis.Visible = 'off';
yticks([0 .1 .2 .3 .4])

% Plot Kir - Capture Speed Low
figure('Position',[500 300 500 500])

% Uncomment for raw traces
for d = 1:length(indlow)
    
    plot((Type(1).Data(indlow(d)).t - Type(1).Data(indlow(d)).t(end)) * 1000,Type(1).Data(indlow(d)).speedVec_Dams(1:length(Type(1).Data(indlow(d)).t)), 'Color',linecolor(1,:));
    hold on
end

shadedErrorBar(((Type(1).Data(3).t - Type(1).Data(3).t(end)) * 1000), meanlVec(1:length(Type(1).Data(3).t)), stdlVec(1:length(Type(1).Data(3).t)),'lineprops',{'Color',linecolor(1,:),'LineWidth',2},'transparent', 1); %shaded area for standard dev


%title('GF1>Kir2.1 - Capture','FontSize',14)
ylabel('Speed (m/s)', 'FontWeight', 'bold','FontSize',14)
xlabel('Time (ms)', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[0,0.325],'Xlim',[-100 0])
set(gca, 'box', 'off')
%ax1 = gca;                   % gca = get current axis
%ax1.YAxis.Visible = 'off';   % remove y-axis
%ax1.XAxis.Visible = 'off';
yticks([0 .1 .2 .3 .4])



% Plot Kir - Capture Speed High & Low combined plot
figure('Position',[500 300 500 500])

% Uncomment for raw traces
for d = 1:length(indhig)
    
    plot((Type(1).Data(indhig(d)).t - Type(1).Data(indhig(d)).t(end)) * 1000,Type(1).Data(indhig(d)).speedVec_Dams(1:length(Type(1).Data(indhig(d)).t)), 'Color',[0, 0.4470, 0.7410],'LineWidth',0.25);
    hold on
end

shadedErrorBar(((Type(1).Data(3).t - Type(1).Data(3).t(end)) * 1000), meanhVec(1:length(Type(1).Data(3).t)), stdhVec(1:length(Type(1).Data(3).t)),'lineprops',{'Color',[0, 0.4470, 0.7410],'LineWidth',1},'transparent', 1); %shaded area for standard dev

hold on

% Uncomment for raw traces
for d = 1:length(indlow)
    
    plot((Type(1).Data(indlow(d)).t - Type(1).Data(indlow(d)).t(end)) * 1000,Type(1).Data(indlow(d)).speedVec_Dams(1:length(Type(1).Data(indlow(d)).t)), 'Color',[0.8500, 0.3250, 0.0980],'LineWidth',0.25);
    hold on
end

shadedErrorBar(((Type(1).Data(3).t - Type(1).Data(3).t(end)) * 1000), meanlVec(1:length(Type(1).Data(3).t)), stdlVec(1:length(Type(1).Data(3).t)),'lineprops',{'Color',[0.8500, 0.3250, 0.0980],'LineWidth',1},'transparent', 1); %shaded area for standard dev



%title('GF1>Kir2.1 - Capture','FontSize',14)
ylabel('Speed (m/s)', 'FontWeight', 'bold','FontSize',14)
xlabel('Time (ms)', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[0,0.325],'Xlim',[-100 0])
set(gca, 'box', 'off')
%ax1 = gca;                   % gca = get current axis
%ax1.YAxis.Visible = 'off';   % remove y-axis
%ax1.XAxis.Visible = 'off';
yticks([0 .1 .2 .3 .4])
%% STEP 31:Individual acceleration plots

%Setup speed variables for plotting
accelp1 = padcat(flipud(Type(1).Data(1).accelVec_Dams), flipud(Type(1).Data(2).accelVec_Dams), flipud(Type(1).Data(3).accelVec_Dams), flipud(Type(1).Data(4).accelVec_Dams), flipud(Type(1).Data(5).accelVec_Dams)...
    ,flipud(Type(1).Data(6).accelVec_Dams), flipud(Type(1).Data(7).accelVec_Dams), flipud(Type(1).Data(8).accelVec_Dams), flipud(Type(1).Data(6).accelVec_Dams), flipud(Type(1).Data(10).accelVec_Dams), flipud(Type(1).Data(11).accelVec_Dams)...
    ,flipud(Type(1).Data(12).accelVec_Dams), flipud(Type(1).Data(13).accelVec_Dams));
accelp1 = flipud(accelp1); %concat accel matricies
numnansp1 = sum(~isnan(accelp1),2); %number of traces averaged as a function of time
meanp1Vec = nanmean(accelp1, 2); %mean accel as a function of time
stdp1Vec = nanstd(accelp1, 0, 2); %standard deviation of accel

accelp2 = padcat(flipud(Type(2).Data(1).accelVec_Dams), flipud(Type(2).Data(2).accelVec_Dams), flipud(Type(2).Data(3).accelVec_Dams), flipud(Type(2).Data(4).accelVec_Dams), flipud(Type(2).Data(5).accelVec_Dams)...
    , flipud(Type(2).Data(6).accelVec_Dams));
accelp2 = flipud(accelp2);
numnansp2 = sum(~isnan(accelp2),2);
meanp2Vec = nanmean(accelp2, 2);
stdp2Vec = nanstd(accelp2, 0, 2);

accelp3 = padcat(flipud(Type(3).Data(1).accelVec_Dams), flipud(Type(3).Data(2).accelVec_Dams), flipud(Type(3).Data(3).accelVec_Dams), flipud(Type(3).Data(4).accelVec_Dams), flipud(Type(3).Data(5).accelVec_Dams));
accelp3 = flipud(accelp3);
numnansp3 = sum(~isnan(accelp3),2);
meanp3Vec = nanmean(accelp3, 2);
stdp3Vec = nanstd(accelp3, 0, 2);

accelp4 = padcat(flipud(Type(4).Data(1).accelVec_Dams), flipud(Type(4).Data(2).accelVec_Dams), flipud(Type(4).Data(3).accelVec_Dams), flipud(Type(4).Data(4).accelVec_Dams), flipud(Type(4).Data(5).accelVec_Dams)...
    ,flipud(Type(4).Data(6).accelVec_Dams), flipud(Type(4).Data(7).accelVec_Dams));
accelp4 = flipud(accelp4);
numnansp4 = sum(~isnan(accelp4),2);
meanp4Vec = nanmean(accelp4, 2);
stdp4Vec = nanstd(accelp4, 0, 2);

linecolor = [[1 0 0] ;[0 0 1] ;[1 .5 0] ;[0 1 0]];

% Plot DL - Capture accel
figure('Position',[500 300 500 500])

% Uncomment for raw traces
for d = 1:length(Type(3).Data)
    
    plot((Type(3).Data(d).t - Type(3).Data(d).t(end)) * 1000,Type(3).Data(d).accelVec_Dams(1:length(Type(3).Data(d).t)), 'Color',linecolor(3,:));
    hold on

    timedata = (Type(3).Data(d).t - Type(3).Data(d).t(end)) * 1000;
    acceldata = Type(3).Data(d).accelVec_Dams(1:length(Type(3).Data(d).t));
    combined = [timedata acceldata];
    a = 0;
end

shadedErrorBar(((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000), meanp3Vec(1:length(Type(3).Data(1).t)), stdp3Vec(1:length(Type(3).Data(1).t)),'lineprops',{'Color',linecolor(3,:),'LineWidth',2},'transparent', 1); %shaded area for standard dev


%title('GF1>+ - Capture','FontSize',14)
ylabel('Acceleration (m/s^{2})', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[-10,7],'Xlim',[-100 0])
yticks([-10 -5 0 5 10])
set(gca, 'box', 'off')
ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';
hold off

% Plot Kir - Capture accel
figure('Position',[500 300 500 500])

% Uncomment for raw traces
for d = 1:length(Type(1).Data)
    
    plot((Type(1).Data(d).t - Type(1).Data(d).t(end)) * 1000,Type(1).Data(d).accelVec_Dams(1:length(Type(1).Data(d).t)), 'Color',linecolor(1,:));
    hold on

    timedata = (Type(1).Data(d).t - Type(1).Data(d).t(end)) * 1000;
    acceldata = Type(1).Data(d).accelVec_Dams(1:length(Type(1).Data(d).t));
    combined = [timedata acceldata];
    a = 0;
end

shadedErrorBar(((Type(1).Data(3).t - Type(1).Data(3).t(end)) * 1000), meanp1Vec(1:length(Type(1).Data(3).t)), stdp1Vec(1:length(Type(1).Data(3).t)),'lineprops',{'Color',linecolor(1,:),'LineWidth',2},'transparent', 1); %shaded area for standard dev


%title('GF1>Kir2.1 - Capture','FontSize',14)
ylabel('Acceleration (m/s^{2})', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[-10,7],'Xlim',[-100 0])
yticks([-10 -5 0 5 10])
set(gca, 'box', 'off')
ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';

% Plot DL - Escape accel
figure('Position',[500 300 500 500])

% Uncomment for raw traces
for d = 1:length(Type(4).Data)
    
    plot((Type(4).Data(d).t - Type(4).Data(d).t(end)) * 1000,Type(4).Data(d).accelVec_Dams(1:length(Type(4).Data(d).t)), 'Color',linecolor(4,:));
    hold on

    timedata = (Type(4).Data(d).t - Type(4).Data(d).t(end)) * 1000;
    acceldata = Type(4).Data(d).accelVec_Dams(1:length(Type(4).Data(d).t));
    combined = [timedata acceldata];
    a = 0;
end

shadedErrorBar(((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000), meanp4Vec(1:length(Type(4).Data(6).t)), stdp4Vec(1:length(Type(4).Data(6).t)),'lineprops',{'Color',linecolor(4,:),'LineWidth',2},'transparent', 1); %shaded area for standard dev


%title('GF1>+ - Escape','FontSize',14)
ylabel('Acceleration (m/s^{2})', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[-10,7],'Xlim',[-100 0])
yticks([-10 -5 0 5 10])
set(gca, 'box', 'off')
ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';
hold off

% Plot Kir - Escape accel
figure('Position',[500 300 500 500])

% Uncomment for raw traces
for d = 1:length(Type(2).Data)
    
    plot((Type(2).Data(d).t - Type(2).Data(d).t(end)) * 1000,Type(2).Data(d).accelVec_Dams(1:length(Type(2).Data(d).t)), 'Color',linecolor(2,:));
    hold on

    timedata = (Type(2).Data(d).t - Type(2).Data(d).t(end)) * 1000;
    acceldata = Type(2).Data(d).accelVec_Dams(1:length(Type(2).Data(d).t));
    combined = [timedata acceldata];
    a = 0;
end

shadedErrorBar(((Type(2).Data(5).t - Type(2).Data(5).t(end)) * 1000), meanp2Vec(1:length(Type(2).Data(5).t)), stdp2Vec(1:length(Type(2).Data(5).t)),'lineprops',{'Color',linecolor(2,:),'LineWidth',2},'transparent', 1); %shaded area for standard dev


%title('GF1>Kir2.1 - Escape','FontSize',14)
ylabel('Acceleration (m/s^{2})', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca,'Ylim',[-10,7],'Xlim',[-100 0])
yticks([-10 -5 0 5 10])
set(gca, 'box', 'off')
ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';


%% STEP 32:Individual n plots

linecolor = [[1 0 0] ;[0 0 1] ;[1 .5 0] ;[0 1 0]];

% Plot DL - Capture Speed
figure('Position',[500 300 500 150])

plot(((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000), numnansp3(1:length(Type(3).Data(1).t)), 'Color',linecolor(3,:)) %number of traces averaged
set(gca,'Ylim', [0, 15], 'Xlim',[-100 0])

xlabel('Time (ms)', 'FontWeight', 'bold','FontSize',14)
ylabel('n', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca, 'box', 'off')
ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';

time1 = ((Type(3).Data(1).t - Type(3).Data(1).t(end)) * 1000);
n1 = numnansp3(1:length(Type(3).Data(1).t));
combined1 = [time1 n1];

% Plot Kir - Capture Speed
figure('Position',[500 300 500 150])

plot(((Type(1).Data(3).t - Type(1).Data(3).t(end)) * 1000), numnansp1(1:length(Type(1).Data(3).t)), 'Color',linecolor(1,:)) %number of traces averaged
set(gca,'Ylim', [0, 15], 'Xlim',[-100 0])

xlabel('Time (ms)', 'FontWeight', 'bold','FontSize',14)
ylabel('n', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca, 'box', 'off')
ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';

time1 = ((Type(1).Data(3).t - Type(1).Data(3).t(end)) * 1000);
n1 = numnansp1(1:length(Type(1).Data(3).t));
combined2 = [time1 n1];

% Plot DL - Escape Speed
figure('Position',[500 300 500 150])

plot(((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000), numnansp4(1:length(Type(4).Data(6).t)), 'Color',linecolor(4,:)) %number of traces averaged
set(gca,'Ylim', [0, 15], 'Xlim',[-100 0])

xlabel('Time (ms)', 'FontWeight', 'bold','FontSize',14)
ylabel('n', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
set(gca, 'box', 'off')
ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';

time1 = ((Type(4).Data(6).t - Type(4).Data(6).t(end)) * 1000);
n1 = numnansp4(1:length(Type(4).Data(6).t));
combined3 = [time1 n1];

% Plot Kir - Escape Speed
figure('Position',[500 300 500 150])

plot(((Type(2).Data(5).t - Type(2).Data(5).t(end)) * 1000), numnansp2(1:length(Type(2).Data(5).t)), 'Color',linecolor(2,:)) %number of traces averaged
set(gca,'Ylim', [0, 15], 'Xlim',[-100 0])

xlabel('Time (ms)', 'FontWeight', 'bold','FontSize',14)
ylabel('n', 'FontWeight', 'bold','FontSize',14)
ax = gca;
ax.FontSize = 14;
ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';
set(gca, 'box', 'off')

time1 = ((Type(2).Data(5).t - Type(2).Data(5).t(end)) * 1000);
n1 = numnansp2(1:length(Type(2).Data(5).t));
combined4 = [time1 n1];

%% STEP 33: Line plots for velocity distribution - based on GF>+ escape
bins = 2;
bin1 = mean(maxSpeed{1,4});
bin2 = max(maxSpeed{1,4});
binind = [0 0.15 Inf];
pctcapDL = [];
pctcapKir = [];

figure('Position',[300 300 500 500])
for i = 1:bins
    for p = 1:4
        index= find(maxSpeed{1,p}>binind(i) & maxSpeed{1,p}<=(binind(i+1)));
        count(p) = length(index);
    end
     pctcapDL(i) = count(4)/(count(3)+count(4));
    pctcapKir(i) = count(2)/(count(1)+count(2));
    countDL(i) = (count(3)+count(4));
    countKir(i) = (count(1)+count(2));
    wsiCapDL_temp = get_error_bars(count(4),(count(3)+count(4)));
    wsiCapKir_temp = get_error_bars(count(2),(count(1)+count(2)));
    wsiCapDL_low(i)=pctcapDL(i)-wsiCapDL_temp(1);
    wsiCapDL_high(i)=wsiCapDL_temp(2)-pctcapDL(i);
    wsiCapKir_low(i)=pctcapKir(i)-wsiCapKir_temp(1);
    wsiCapKir_high(i)=wsiCapKir_temp(2)-pctcapKir(i);
end
x = 1:bins;
b = plot(x,100*[pctcapDL' pctcapKir'],'LineWidth',1.5);
b(1).Color = [0 0 0];
b(2).Color = [.75 .75 .75];
b(1).Marker = 'o';
b(2).Marker = 'o';
b(1).MarkerFaceColor = [0 0 0];
b(2).MarkerFaceColor = [.75 .75 .75];

ylabel('Percent Escape (%)', 'FontWeight', 'bold')
xlabel('Peak Speed (m/s)', 'FontWeight', 'bold')
legend('GF1>+','GF1>Kir2.1','Location', 'northeast')
 ylim([0 100])
% xlim([2.5 4.5])
ax = gca;
ax.FontSize = 14;

xticks(x)
xticklabels({'<0.15','>0.15'})
set(gca, 'box', 'off')

figure('Position',[300 300 500 500])
b = bar(x,100*[pctcapDL' pctcapKir'],'LineWidth',1.5);
b(1).FaceColor = [.5 .5 .5];
b(2).FaceColor = [1 1 1];
hold on
er = errorbar([.86,1.14;1.86,2.14],100*[pctcapDL' pctcapKir'],100*[wsiCapDL_low' wsiCapKir_low'],100*[wsiCapDL_high' wsiCapKir_high'],'LineStyle','none','Color',[0 0 0],'LineWidth',1.5);    
%er.Color = [0 0 0];                            
%er.LineStyle = 'none'; 
%errorbar(wsiCapDL,wsiCapKir)

ylabel('Percent Escape (%)', 'FontWeight', 'bold')
xlabel('Peak Speed (m/s)', 'FontWeight', 'bold')
legend('GF1>+','GF1>Kir2.1','Location', 'northeast')
 ylim([0 100])
% xlim([2.5 4.5])
ax = gca;
ax.FontSize = 14;

xticks(x)
xticklabels({'<0.15','>0.15'})
set(gca, 'box', 'off')

%% STEP 34: Max Speed Histograms and Bar Plots

dlmaxspeed=[maxSpeed{3};maxSpeed{4}];
kirmaxspeed=[maxSpeed{1};maxSpeed{2}];

dlavgspeed=[avgSpeed{3};avgSpeed{4}];
kiravgspeed=[avgSpeed{1};avgSpeed{2}];

dlmaxacc=[peakAcc{3};peakAcc{4}];
kirmaxacc=[peakAcc{1};peakAcc{2}];


dltimeorder = [6 1 7 8 9 2 10 3 11 4 12 5];
kirtimeorder = [1 14 2 3 4 15 16 5 6 17 7 8 9 10 11 12 13 18 19];

dltimestamp = [0 87 177 204 244 281 304 314 331 414 422 447]; % first event is t=0, time elapsed in min
kirtimestamp= [0 13 13 13 26 75 75 83 88 88 100 217 258 339 355 367 428 466 476];

dlspeedtime = dlmaxspeed(dltimeorder); %ordering max speeds by time video was captured
kirspeedtime = kirmaxspeed(kirtimeorder);

dlaspeedtime = dlavgspeed(dltimeorder); %ordering avg speeds by time video was captured
kiraspeedtime = kiravgspeed(kirtimeorder);

dlacctime = dlmaxacc(dltimeorder); %ordering max acc by time video was captured
kiracctime = kirmaxacc(kirtimeorder);

numbins = 4;
bins=min(dlmaxspeed):(max(dlmaxspeed)-min(dlmaxspeed))/numbins:max(dlmaxspeed);
figure
histogram(dlmaxspeed,bins) %,'Normalization','Probability')
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
xlabel('GF1>+ Damselfly Max Speed')
ylabel('Percentage Damselflies')
%ylim([0,7])
axisHandle = gca;
histHandle = axisHandle.Children;
histdl = histHandle.Values; 

figure
histogram(kirmaxspeed,bins)%,'Normalization','Probability')
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
xlabel('GF1>Kir2.1 Damselfly Max Speed')
ylabel('Percentage Damselflies')
%ylim([0,7])
axisHandle = gca;
histHandle = axisHandle.Children;
histkir = histHandle.Values; 

binmids = (bins(1:end-1)+bins(2:end))./2;

figure('Position',[300 300 500 500])
b = bar(binmids,100*[histdl'/length(dlmaxspeed)  histkir'/length(kirmaxspeed)],'LineWidth',1.5);
b(1).FaceColor = [.5 .5 .5];
b(2).FaceColor = [1 1 1];

ylabel('Percentage Damselflies', 'FontWeight', 'bold')
xlabel('Peak Speed (m/s)', 'FontWeight', 'bold')
legend('GF1>+','GF1>Kir2.1','Location', 'northeast')
 ylim([0 40])
% xlim([2.5 4.5])
ax = gca;
ax.FontSize = 14;

%xticks(x)
%xticklabels({'<0.13','>0.13'})
set(gca, 'box', 'off')

b = plot(binmids,100*[histdl'/length(dlmaxspeed)  histkir'/length(kirmaxspeed)],'LineWidth',1.5);
b(1).Color = [0 0 0];
b(2).Color = [.75 .75 .75];
b(1).Marker = 'o';
b(2).Marker = 'o';
b(1).MarkerFaceColor = [0 0 0];
b(2).MarkerFaceColor = [.75 .75 .75];

ylabel('Percentage Damselflies','FontWeight','bold')
xlabel('Peak Speed (m/s)', 'FontWeight', 'bold')
legend('GF1>+','GF1>Kir2.1','Location', 'northeast')
 ylim([10 40])
% xlim([2.5 4.5])
ax = gca;
ax.FontSize = 14;

%xticks(x)
%xticklabels({'<0.24','>0.24'})
set(gca, 'box', 'off')


figure
xdata = repmat(1, length(dlmaxspeed), 1);
scatter(xdata, dlmaxspeed, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7]);
hold on
plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(dlmaxspeed, 1), 2, 1), 'k-','LineWidth',3)
hold on
plot(repmat(xdata(1),2,1),[mean(dlmaxspeed, 1)-(std(dlmaxspeed)/sqrt(length(dlmaxspeed)));mean(dlmaxspeed, 1)+(std(dlmaxspeed)/sqrt(length(dlmaxspeed)))],'k-','LineWidth',1.5)
hold on

xdata = repmat(2, length(kirmaxspeed), 1);
scatter(xdata, kirmaxspeed, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7]);
hold on
plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(kirmaxspeed, 1), 2, 1), 'k-','LineWidth',3)
hold on
plot(repmat(xdata(1),2,1),[mean(kirmaxspeed, 1)-(std(kirmaxspeed)/sqrt(length(kirmaxspeed)));mean(kirmaxspeed, 1)+(std(kirmaxspeed)/sqrt(length(kirmaxspeed)))],'k-','LineWidth',1.5)
hold on

ylabel('Peak Speed (m/s)', 'FontWeight', 'bold')
ylim([0 0.35])
xlim([.5 2.5])
ax = gca;
ax.FontSize = 14;

xticks([1 2])
yticks([0 .1 .2 .3 .4])
xticklabels({'GF1>+','GF1>Kir2.1'})

%% Time ordered metrics

figure
plot(dltimestamp,dlspeedtime,'ko','MarkerFaceColor','k')
xlabel('Time Elapsed (min)')
ylabel('Peak Speed (m/s)')
title('GF>+')
ylim([0,0.35])
xlim([0 480])

figure
plot(kirtimestamp,kirspeedtime,'ko','MarkerFaceColor','k')
xlabel('Time Elapsed (min)')
ylabel('Peak Speed (m/s)')
title('GF>Kir2.1')
ylim([0,0.35])
xlim([0 480])

figure
plot(dltimestamp,dlaspeedtime,'ko','MarkerFaceColor','k')
xlabel('Time Elapsed (min)')
ylabel('Average Speed (m/s)')
title('GF>+')
ylim([0,0.25])
xlim([0 480])

figure
plot(kirtimestamp,kiraspeedtime,'ko','MarkerFaceColor','k')
xlabel('Time Elapsed (min)')
ylabel('Average Speed (m/s)')
title('GF>Kir2.1')
ylim([0,0.25])
xlim([0 480])

figure
plot(dltimestamp,dlacctime,'ko','MarkerFaceColor','k')
xlabel('Time Elapsed (min)')
ylabel('Peak Acceleration (m/s^{2})')
title('GF>+')
ylim([0,6])
xlim([0 480])

figure
plot(kirtimestamp,kiracctime,'ko','MarkerFaceColor','k')
xlabel('Time Elapsed (min)')
ylabel('Peak Acceleration (m/s^{2})')
title('GF>Kir2.1')
ylim([0,6])
xlim([0 480])

%% Attack Duration

linecolor = [[1 0 0] ;[0 0 1] ;[1 .5 0] ;[0 1 0]];

figure('Position',[300 300 600 500])
for p = 1:4
    
    attackDur = zeros(length(Type(p).Data),1);
    
    for j = 1:length(Type(p).Data)
    attackDur(j) = length(Type(p).Data(j).t);
    end
    
    index = [2 4 1 3];
    xdata = repmat(index(p), length(attackDur), 1);
    scatter(xdata, attackDur, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',linecolor(p,:),'MarkerEdgeColor',linecolor(p,:));
    hold on
    
    plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(attackDur), 2, 1), 'k-','LineWidth',3)
    hold on
    plot(repmat(xdata(1),2,1),[mean(attackDur)-(std(attackDur)/sqrt(length(attackDur)));mean(attackDur)+(std(attackDur)/sqrt(length(attackDur)))],'k-','LineWidth',1.5)
    hold on
end

%title('GF1>+','FontSize',12)
ylabel('Attack Duration (ms)', 'FontWeight', 'bold')
% xlabel('Fly Genotype', 'FontWeight', 'bold')
%ylim([-40 40])
xlim([.5 4.5])
ax = gca;
ax.FontSize = 14;

xticks([1 2 3 4])
xticklabels({'GF1>+','GF1>Kir2.1','GF1>+','GF1>Kir2.1'})
xlabel('Capture                             Escape')

%% Angular Size at Start of Attack

linecolor = [[1 0 0] ;[0 0 1] ;[1 .5 0] ;[0 1 0]];

figure('Position',[300 300 600 500])
for p = 1:4
    
    angSizeStart = zeros(length(Type(p).Data),1);
    
    for j = 1:length(Type(p).Data)
    angSizeStart(j) = Type(p).Data(j).angSizetable(1);
    end
    
    index = [2 4 1 3];
    xdata = repmat(index(p), length(angSizeStart), 1);
    scatter(xdata, angSizeStart, 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',linecolor(p,:),'MarkerEdgeColor',linecolor(p,:));
    hold on
    
    plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(angSizeStart), 2, 1), 'k-','LineWidth',3)
    hold on
    plot(repmat(xdata(1),2,1),[mean(angSizeStart)-(std(angSizeStart)/sqrt(length(angSizeStart)));mean(angSizeStart)+(std(angSizeStart)/sqrt(length(angSizeStart)))],'k-','LineWidth',1.5)
    hold on
end

%title('GF1>+','FontSize',12)
ylabel('Starting Angular Size (°)', 'FontWeight', 'bold')
% xlabel('Fly Genotype', 'FontWeight', 'bold')
%ylim([-40 40])
xlim([.5 4.5])
ax = gca;
ax.FontSize = 14;

xticks([1 2 3 4])
xticklabels({'GF1>+','GF1>Kir2.1','GF1>+','GF1>Kir2.1'})
xlabel('Capture                             Escape')

%% Starting distance vs speed/accel

linecolor = [[1 0 0] ;[0 0 1] ;[1 .5 0] ;[0 1 0]];

figure('Position',[300 300 600 500])
for p = 1:4
    index = [3 1 4 2];
    distStart = zeros(length(Type(index(p)).Data),1);
    
    for j = 1:length(Type(index(p)).Data)
    distStart(j) = 1000*Type(index(p)).Data(j).distVec_Dams(8);
    end
    
    scatter(distStart', maxSpeed{index(p)},'MarkerFaceColor',linecolor(index(p),:),'MarkerEdgeColor',linecolor(index(p),:));
    hold on
end

%title('GF1>+','FontSize',12)
ylabel('Peak Speed (m/s)', 'FontWeight', 'bold')
xlabel('Starting Distance (mm)', 'FontWeight', 'bold')
ylim([0 .35])
%xlim([.5 4.5])
ax = gca;
ax.FontSize = 14;
legend('GF1>+ - Capture','GF1>Kir2.1 - Capture','GF1>+ - Escape','GF1>Kir2.1 - Escape','Location', 'northeast')

figure('Position',[300 300 600 500])
for p = 1:4
    index = [3 1 4 2];
    distStart = zeros(length(Type(index(p)).Data),1);
    
    for j = 1:length(Type(index(p)).Data)
    distStart(j) = 1000*Type(index(p)).Data(j).distVec_Dams(8);
    end
    
    scatter(distStart', peakAcc{index(p)},'MarkerFaceColor',linecolor(index(p),:),'MarkerEdgeColor',linecolor(index(p),:));
    hold on
end

%title('GF1>+','FontSize',12)
ylabel('Peak Acceleration (m/s^{2})', 'FontWeight', 'bold')
xlabel('Starting Distance (mm)', 'FontWeight', 'bold')
ylim([0 6])
%xlim([.5 4.5])
ax = gca;
ax.FontSize = 14;
legend('GF1>+ - Capture','GF1>Kir2.1 - Capture','GF1>+ - Escape','GF1>Kir2.1 - Escape','Location', 'northeast')

%% Attack duration vs speed/accel/distance

linecolor = [[1 0 0] ;[0 0 1] ;[1 .5 0] ;[0 1 0]];

figure('Position',[300 300 600 500])
for p = 1:4
    index = [3 1 4 2];
    attackDur = zeros(length(Type(index(p)).Data),1);
    
    for j = 1:length(Type(index(p)).Data)
        attackDur(j) = length(Type(index(p)).Data(j).t);
    end
    
    scatter(maxSpeed{index(p)},attackDur,'MarkerFaceColor',linecolor(index(p),:),'MarkerEdgeColor',linecolor(index(p),:));
    hold on
end

%title('GF1>+','FontSize',12)
xlabel('Peak Speed (m/s)', 'FontWeight', 'bold')
ylabel('Attack Duration (ms)', 'FontWeight', 'bold')
%ylim([0 .35])
%xlim([.5 4.5])
ax = gca;
ax.FontSize = 14;
%legend('GF1>+ - Capture','GF1>Kir2.1 - Capture','GF1>+ - Escape','GF1>Kir2.1 - Escape','Location', 'northeast')

figure('Position',[300 300 600 500])
for p = 1:4
    index = [3 1 4 2];
    attackDur = zeros(length(Type(index(p)).Data),1);
    
    for j = 1:length(Type(index(p)).Data)
        attackDur(j) = length(Type(index(p)).Data(j).t);
    end
    
    scatter(peakAcc{index(p)},attackDur,'MarkerFaceColor',linecolor(index(p),:),'MarkerEdgeColor',linecolor(index(p),:));
    hold on
end

%title('GF1>+','FontSize',12)
xlabel('Peak Acceleration (m/s^{2})', 'FontWeight', 'bold')
ylabel('Attack Duration (ms)', 'FontWeight', 'bold')
%ylim([0 6])
%xlim([.5 4.5])
ax = gca;
ax.FontSize = 14;
%legend('GF1>+ - Capture','GF1>Kir2.1 - Capture','GF1>+ - Escape','GF1>Kir2.1 - Escape','Location', 'northeast')

figure('Position',[300 300 600 500])
for p = 1:4
    index = [3 1 4 2];
    distStart = zeros(length(Type(index(p)).Data),1);
    
    for j = 1:length(Type(index(p)).Data)
    distStart(j) = 1000*Type(index(p)).Data(j).distVec_Dams(8);
    end
    
    attackDur = zeros(length(Type(index(p)).Data),1);
    
    for j = 1:length(Type(index(p)).Data)
        attackDur(j) = length(Type(index(p)).Data(j).t);
    end
    
    scatter(distStart',attackDur,'MarkerFaceColor',linecolor(index(p),:),'MarkerEdgeColor',linecolor(index(p),:));
    hold on
end

%title('GF1>+','FontSize',12)
ylabel('Attack Duration (ms)', 'FontWeight', 'bold')
xlabel('Starting Distance (mm)', 'FontWeight', 'bold')
%ylim([0 6])
%xlim([.5 4.5])
ax = gca;
ax.FontSize = 14;
%legend('GF1>+ - Capture','GF1>Kir2.1 - Capture','GF1>+ - Escape','GF1>Kir2.1 - Escape','Location', 'northeast')


%% Loom Angular Size over Time

initStimSize = 10;
finalStimSize = 180;
ellovervee = 10;
stimTimeStep = (1/360)*1000;%milliseconds per frame channel at 120 Hz

%making disk radius per time (frame) vector
minTheta = deg2rad(initStimSize);
maxTheta = deg2rad(finalStimSize);
stimStartTime = ellovervee/tan(minTheta/2);
stimEndTime = ellovervee/tan(maxTheta/2);
stimTotalDuration = stimStartTime-stimEndTime;
stimTimeVector = fliplr(stimEndTime:stimTimeStep:stimStartTime);
stimThetaVector = 2.*atan(ellovervee./stimTimeVector);
stimThetaVector = rad2deg(stimThetaVector);

time = 0:stimTotalDuration/(length(stimThetaVector)-1):stimTotalDuration;

figure('Position',[1000,500,500,250])
plot(time,stimThetaVector,'-k','LineWidth',2)
hold on
plot([0,0],[0,stimThetaVector(1)],'-k','LineWidth',2)
hold on 
plot([stimTotalDuration,stimTotalDuration+100],[stimThetaVector(end),stimThetaVector(end)],'-k','LineWidth',2)
xlabel('Time (ms)')
ylabel('Loom Angular Size (deg)')
set(gca, 'box', 'off')
xlim([0,500])

simvector = [time stimThetaVector];

%% Wing Lift to Takeoff Ethograms

%Wing Lift Times
GF110DLWL = rmmissing((cell2mat(GF1DL10.frame_of_wing_movement) - GF1DL10.stimStart) /6);
GF140DLWL = rmmissing((cell2mat(GF1DL40.frame_of_wing_movement) - GF1DL40.stimStart) /6);
GF210DLWL = rmmissing((cell2mat(GF2DL10.frame_of_wing_movement) - GF2DL10.stimStart) /6);
GF240DLWL = rmmissing((cell2mat(GF2DL40.frame_of_wing_movement) - GF2DL40.stimStart) /6);
DL10WL = [GF110DLWL; GF210DLWL];
DL40WL = [GF140DLWL; GF240DLWL];

indDL10 = find(DL10WL<=0);
DL10WL(indDL10) = [];
indDL40 = find(DL40WL<=0);
DL40WL(indDL40) = [];

GF110KIRWL = rmmissing((cell2mat(GF1Kir10.frame_of_wing_movement) - GF1Kir10.stimStart) /6);
GF140KIRWL = rmmissing((cell2mat(GF1Kir40.frame_of_wing_movement) - GF1Kir40.stimStart) /6);
GF210KIRWL = rmmissing((cell2mat(GF2Kir10.frame_of_wing_movement) - GF2Kir10.stimStart) /6);
GF240KIRWL = rmmissing((cell2mat(GF2Kir40.frame_of_wing_movement) - GF2Kir40.stimStart) /6);
KIR10WL = [GF110KIRWL; GF210KIRWL];
KIR40WL = [GF140KIRWL; GF240KIRWL];

indKIR10 = find(KIR10WL<=0);
KIR10WL(indKIR10) = [];
indKIR40 = find(KIR40WL<=0);
KIR40WL(indKIR40) = [];


%Take off Times
GF110DLTO = rmmissing((cell2mat(GF1DL10.frame_of_take_off) - GF1DL10.stimStart) /6);
GF140DLTO = rmmissing((cell2mat(GF1DL40.frame_of_take_off) - GF1DL40.stimStart) /6);
GF210DLTO = rmmissing((cell2mat(GF2DL10.frame_of_take_off) - GF2DL10.stimStart) /6);
GF240DLTO = rmmissing((cell2mat(GF2DL40.frame_of_take_off) - GF2DL40.stimStart) /6);
DL10TO = [GF110DLTO; GF210DLTO];
DL40TO = [GF140DLTO; GF240DLTO];
% 
 DL10TO(indDL10) = [];
 DL40TO(indDL40) = [];

GF110KIRTO = rmmissing((cell2mat(GF1Kir10.frame_of_take_off) - GF1Kir10.stimStart) /6);
GF140KIRTO = rmmissing((cell2mat(GF1Kir40.frame_of_take_off) - GF1Kir40.stimStart) /6);
GF210KIRTO = rmmissing((cell2mat(GF2Kir10.frame_of_take_off) - GF2Kir10.stimStart) /6);
GF240KIRTO = rmmissing((cell2mat(GF2Kir40.frame_of_take_off) - GF2Kir40.stimStart) /6);
KIR10TO = [GF110KIRTO; GF210KIRTO];
KIR40TO = [GF140KIRTO; GF240KIRTO];

 KIR10TO(indKIR10) = [];
 KIR40TO(indKIR40) = [];

[DL10WL,index]=sort(DL10WL,'descend');
DL10TO=DL10TO(index);
[DL40WL,index]=sort(DL40WL,'descend');
DL40TO=DL40TO(index);

[KIR10WL,index]=sort(KIR10WL,'descend');
KIR10TO=KIR10TO(index);
[KIR40WL,index]=sort(KIR40WL,'descend');
KIR40TO=KIR40TO(index);

%Reduce trials to match number for 1/v 10 GF>Kir trials
for i = 1:(length(DL10WL)-length(KIR10WL))
    check = 1;
    while check == 1
        ind = round(rand*length(DL10WL));
        if ind == 0
            ind = 1;
        end
        if ~isnan(DL10WL(ind))
            DL10WL(ind) = NaN;
            check = 0;
        end
    end
end

ind2 = find(isnan(DL10WL));
DL10WL(ind2) = [];
DL10TO(ind2) = [];

for i = 1:(length(DL40WL)-length(KIR10WL))
    check = 1;
    while check == 1
        ind = round(rand*length(DL40WL));
        if ind == 0
            ind = 1;
        end
        if ~isnan(DL40WL(ind))
            DL40WL(ind) = NaN;
            check = 0;
        end
    end
end

ind2 = find(isnan(DL40WL));
DL40WL(ind2) = [];
DL40TO(ind2) = [];

for i = 1:(length(KIR40WL)-length(KIR10WL))
    check = 1;
    while check == 1
        ind = round(rand*length(KIR40WL));
        if ind == 0
            ind = 1;
        end
        if ~isnan(KIR40WL(ind))
            KIR40WL(ind) = NaN;
            check = 0;
        end
    end
end

ind2 = find(isnan(KIR40WL));
KIR40WL(ind2) = [];
KIR40TO(ind2) = [];


figure('Position',[1000,500,500,300])
for i = 1:length(DL10WL)
    plot ([DL10WL(i),DL10TO(i)],[i,i],'k','LineWidth',1)
    hold on
end
xlim([0,150])
ylim([0,30])
xlabel('time (ms)')
ylabel('GF>+ Trials')
set(gca, 'box', 'off')

figure('Position',[1000,500,500,300])
for i = 1:length(DL40WL)
    plot ([DL40WL(i),DL40TO(i)],[i,i],'k','LineWidth',1)
    hold on
end
xlim([0,500])
ylim([0,30])
xlabel('time (ms)')
ylabel('GF>+ Trials')
set(gca, 'box', 'off')

figure('Position',[1000,500,500,300])
for i = 1:length(KIR10WL)
    plot ([KIR10WL(i),KIR10TO(i)],[i,i],'k','LineWidth',1)
    hold on
end
xlim([0,150])
ylim([0,30])
xlabel('time (ms)')
ylabel('GF>Kir2.1 Trials')
set(gca, 'box', 'off')

figure('Position',[1000,500,500,300])
for i = 1:length(KIR40WL)
    plot ([KIR40WL(i),KIR40TO(i)],[i,i],'k','LineWidth',1)
    hold on
end
xlim([0,500])
ylim([0,30])
xlabel('time (ms)')
ylabel('GF>Kir2.1 Trials')
set(gca, 'box', 'off')

%% Box plots for start of attack to capture/escape
linecolor = [[1 0 0] ;[0 0 1] ;[1 .5 0] ;[0 1 0]];


figure('Position',[300 300 500 500])
for p = 1:4
    index = [2 4 1 3];
    xdata = repmat(index(p), length(eventdur{1,p}'), 1);
    scatter(xdata, eventdur{1,p}', 'o', 'jitter','on', 'jitterAmount', 0.05,'MarkerFaceColor',linecolor(p,:),'MarkerEdgeColor',linecolor(p,:));
    hold on
    
    plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(eventdur{1,p}', 1), 2, 1), 'k-','LineWidth',3)
    hold on
    plot(repmat(xdata(1),2,1),[mean(eventdur{1,p}', 1)-(std(eventdur{1,p}')/sqrt(length(eventdur{1,p}')));mean(eventdur{1,p}', 1)+(std(eventdur{1,p}')/sqrt(length(eventdur{1,p}')))],'k-','LineWidth',1.5)
    hold on
end




%title('GF1 > +','FontSize',12)
ylabel('Damselfly Attack Duration (ms)', 'FontWeight', 'bold')
%xlabel('Fly Genotype', 'FontWeight', 'bold')
%ylim([0 120])
xlim([.5 4.5])
ax = gca;
ax.FontSize = 14;

% labels = {'line1 line2','line1 line2','line1 line2'};
% labels = cellfun(@(x) strrep(x,' ','\newline'), labels,'UniformOutput',false);
% a = gca;
% a.XTickLabel = labels;


xticks([1 2 3 4])
xticklabels({'GF1>+','GF1>Kir2.1','GF1>+','GF1>Kir2.1'})
xlabel('Capture                      Escape')
% set(gca,'Ydir','reverse')
% set(gca,'yscale','log')


% a2 = axes('YAxisLocation', 'Right');
% % Hide second plot.
% set(a2, 'color', 'none')
% set(a2, 'XTick', [])
% % Set scala for second Y.
% set(a2, 'YLim', [20 25])
