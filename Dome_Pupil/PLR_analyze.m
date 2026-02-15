%% Pupillary Light Reflex (PLR) Analysis - Dome Setup
% This script analyzes pupillary responses to different illuminance levels
% in a dome/projection setup. It fits sigmoidal dose-response curves to
% extract EC50 (half-maximal effective concentration) and Hill slope values.
%
% Input files (selected via GUI):
%   - Parameter file: Contains stimIn/stimOut structures with experiment parameters
%   - Left eye trace: Pupil area measurements over time
%   - Right eye trace: Pupil area measurements over time
%
% Outputs:
%   - Figures: Dose-response curves with sigmoidal fits
%   - Saved .mat file: All analysis results including EC50, Hill slope, fitted curves
%
% Author: Fitzpatrick et al., 2024

clear; clc; close all;

%% Define variables and load data

% Experimental parameters - adjust these for your setup
stim_fps = 60;  % Stimulus presentation frame rate (Hz)
movie_fps = 15; % Pupil tracking camera frame rate (Hz)
mm_per_pix = 0.00345; % Conversion: millimeters per pixel for pupil area measurements
                      % (depends on camera and magnification)

% Illuminance conversion
conv = 1; % Set to 1 to convert cd/m² to R* (photoisomerizations)
cd_to_R = 13.8893; % Empirical conversion factor between cd/m² and R*
                   % (photoisomerizations per photoreceptor per second)

%parameters
[fileName, pathName]=uigetfile ('*.mat','Select parameter data');
load([pathName fileName]);

if conv==1
    
    stimIn.ExpLum = stimIn.ExpLum.*cd_to_R;
    
end

for i = 1:length(stimOut)
    for j = 1:6
        Data{i,j} = strtrim(convertCharsToStrings(stimOut(i,j,:)));
        if j ~= 3
            Data{i,j} = str2double(Data{i,j});
        else
            if Data{i,j}=="Adapt"
                Data{i,j}=0;
            elseif Data{i,j} == "Pre"
                Data{i,j} = 0.5;
            elseif Data{i,j} == "Stim"
                Data{i,j} = 1;
            elseif Data{i,j} == "Post"
                Data{i,j} = 0.5;
            end
        end
    end
end

Data = cell2mat(Data);

disp(['Succesfully loaded ' fileName])

%traces
[lfileName, lpathName]=uigetfile ('*.mat','Select Left Eye Trace');
if lfileName ~=0
    load([lpathName lfileName]);
    disp(['Succesfully loaded ' lfileName])
    LArea = area(1:length(Data)/stim_fps*movie_fps)*(mm_per_pix)^2;
else
    LArea = NaN(length(Data)/stim_fps*movie_fps,1); %in case there's no file
end

[rfileName, lpathName]=uigetfile ('*.mat','Select Right Eye Trace');
if rfileName ~=0
    load([lpathName rfileName]);
    disp(['Succesfully loaded ' rfileName])
    RArea = area(1:length(Data)/stim_fps*movie_fps)*(mm_per_pix)^2;
else
    RArea = NaN(length(Data)/stim_fps*movie_fps,1); %in case there's no file
end

%% speed

% Extract speed/distance of mouse from wheel encoder
% Wheel parameters: 15 cm diameter, 4028 pulses per rotation
mSpeed = zeros(length(corrData),1);
mSpeed(2:end) = abs(diff(corrData(:,end)));
mSpeed(mSpeed>150 | mSpeed<-150)= 0; % Remove points where encoder cycles (Arduino bit limit)
% Convert to cm/s: (pulses/frame) * (circumference/pulses_per_rotation) * fps
mSpeed = mSpeed*(15*pi*stim_fps/4028); % now in cm/s
mSpeed = downsample(mSpeed,stim_fps/movie_fps); % Match pupil tracking frame rate
smoothSpeed = movmean(mSpeed,movie_fps/2,'omitnan'); % Smooth with 0.5s window

%% Process data
time = (tEnd-tStart)/60;
Data(:,1) = linspace(1/stim_fps,time,length(Data));
runTime = double(stimIn.Adapt+stimIn.Pre+stimIn.Dur+stimIn.Post)*stim_fps;
PupTime = double(runTime/stim_fps*movie_fps);

LArea_runs = reshape(LArea,PupTime,length(stimIn.ExpLum));
LArea_pre = mean(LArea_runs(stimIn.Adapt*movie_fps:(stimIn.Adapt+stimIn.Pre)*movie_fps,:),'omitnan');

RArea_runs = reshape(RArea,PupTime,length(stimIn.ExpLum));
RArea_pre = mean(RArea_runs(stimIn.Adapt*movie_fps:(stimIn.Adapt+stimIn.Pre)*movie_fps,:),'omitnan');

smoothLPupil = movmean(LArea_runs,movie_fps/2,'omitnan');
[L_M, L_I] = min(smoothLPupil(((stimIn.Adapt+stimIn.Pre+stimIn.Dur/2)*movie_fps:(stimIn.Adapt+stimIn.Pre+stimIn.Dur)*movie_fps),:));
L_I = L_I+double((stimIn.Adapt+stimIn.Pre+stimIn.Dur/2)*movie_fps);

smoothRPupil = movmean(RArea_runs,movie_fps/2,'omitnan');
[R_M, R_I] = min(smoothRPupil(((stimIn.Adapt+stimIn.Pre+stimIn.Dur/2)*movie_fps:(stimIn.Adapt+stimIn.Pre+stimIn.Dur)*movie_fps),:));
R_I = R_I+double((stimIn.Adapt+stimIn.Pre+stimIn.Dur/2)*movie_fps);

%% Calculate

% Calculate pupil area at minimum constriction for each trial
PupilChangeLeft = zeros(1,length(stimIn.ExpLum));
PupilChangeRight = zeros(1,length(stimIn.ExpLum));

for i = 1:length(stimIn.ExpLum)
   
    % Average pupil area in ±2.5s window around minimum
    if isnan(L_I(i))
        PupilChangeLeft(i) = NaN;
    else
        PupilChangeLeft(i) = mean(smoothLPupil((L_I(i)-round(2.5*movie_fps)):(L_I(i)+round(2.5*movie_fps)),i),'omitnan');
    end
    
    if isnan(R_I(i))
        PupilChangeRight(i) = NaN;
    else
        PupilChangeRight(i) = mean(smoothRPupil((R_I(i)-round(2.5*movie_fps)):(R_I(i)+round(2.5*movie_fps)),i),'omitnan');
    end
end

%% Plot

% Sigmoidal function for dose-response curve fitting
% Parameters: b(1)=minimum, b(2)=maximum, b(3)=EC50, b(4)=Hill slope
% Note: EC50 is the illuminance at halfway between min and max response
SigmoidFit = @(b,x)(b(1)+(b(2)-b(1))./(1+(b(3)./x).^b(4)));
%Initial guess for least squares regression
beta0 = [0,0.5,10,-1];
%define lower and upper bounds. 
lb = [0,0,0.00001,-10];
ub = [0.1,1,10000000000,0];


%Plot and Fit for Right
figure
points = semilogx(stimIn.ExpLum,PupilChangeRight, 'ko');
points.MarkerSize = 7;
points.MarkerFaceColor = [0 0 0];
ax = gca;
ax.Box = 'off';
ax.YLabel.String = 'Pupil Area (mm^2)';
ax.YLabel.FontSize = 12;
ax.YLim = [0 0.7];
ax.YTick = [0.0 0.1 0.2 0.3 0.4 0.5 0.6];
ax.XLim = [stimIn.ExpLum(1)/10 stimIn.ExpLum(end)*10];
ax.XTick = [10^-1 10^0 10^1 10^2 10^3];
ax.XLabel.String = 'Illuminance (Rod cd/m^2)';
FitPupilChangeR = lsqcurvefit(SigmoidFit,beta0,stimIn.ExpLum,PupilChangeRight,lb,ub);
R_EC50 = FitPupilChangeR(3);
R_Hill = FitPupilChangeR(4);
hold on
logrange = logspace(log10(stimIn.ExpLum(1)/10), log10(stimIn.ExpLum(end)*10));
line = semilogx(logrange, SigmoidFit(FitPupilChangeR,logrange));
line.Color = 'k';
line.LineWidth = 1;
maxref = refline(0,FitPupilChangeR(2));
maxref.Color = 'k';
maxref.LineWidth = 1;
maxref.LineStyle = '--';
minref = refline(0,FitPupilChangeR(1));
minref.Color = 'k';
minref.LineWidth = 1;
minref.LineStyle = '--';

%Plot and Fit for Left
figure
points = semilogx(stimIn.ExpLum,PupilChangeLeft, 'ko');
points.MarkerSize = 7;
points.MarkerFaceColor = [0 0 0];
ax = gca;
ax.Box = 'off';
ax.YLabel.String = 'Pupil Area (mm^2)';
ax.YLabel.FontSize = 12;
ax.YLim = [0 0.7];
ax.YTick = [0.0 0.1 0.2 0.3 0.4 0.5 0.6];
ax.XLim = [stimIn.ExpLum(1)/10 stimIn.ExpLum(end)*10];
ax.XTick = [10^-1 10^0 10^1 10^2 10^3];
ax.XLabel.String = 'Illuminance (Rod cd/m^2)';
FitPupilChangeL = lsqcurvefit(SigmoidFit,beta0,stimIn.ExpLum,PupilChangeLeft,lb,ub);
L_EC50 = FitPupilChangeL(3);
L_Hill = FitPupilChangeL(4);
hold on
logrange = logspace(log10(stimIn.ExpLum(1)/10), log10(stimIn.ExpLum(end)*10));
line = semilogx(logrange, SigmoidFit(FitPupilChangeL,logrange));
line.Color = 'k';
line.LineWidth = 1;
maxref = refline(0,FitPupilChangeL(2));
maxref.Color = 'k';
maxref.LineWidth = 1;
maxref.LineStyle = '--';
minref = refline(0,FitPupilChangeL(1));
minref.Color = 'k';
minref.LineWidth = 1;
minref.LineStyle = '--';

%% Save
[save_fileName, save_pathName]=uiputfile ('*.mat','Save PLR Analysis');
save([save_pathName save_fileName]);

