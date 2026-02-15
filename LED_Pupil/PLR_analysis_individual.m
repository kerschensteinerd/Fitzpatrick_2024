%% Pupillary Light Reflex (PLR) Analysis - LED Setup (Individual)
% This script analyzes pupillary responses to LED stimuli for individual subjects.
% It fits sigmoidal dose-response curves to extract EC50 and Hill slope values.
%
% Input files (selected via GUI):
%   - Parsed data file: Contains pupil traces and stimulus parameters
%
% Outputs:
%   - Figure: Dose-response curve with sigmoidal fit
%   - Saved .mat file: EC50, Hill slope, and all analysis results
%
% Author: Fitzpatrick et al., 2024

clear; clc; close all;

%%Load Parsed files
[fileName, pathName] = uigetfile('.mat', 'Load parsed file');
%Load first file to get parameters and size array
load([pathName fileName]);
disp(['Succesfully loaded ' fileName])

short = 0; %=1 if short and does not include post stim time;

%resort by RStar in case stimulus program messed this up
sortedAllCalib = sortrows(AllCalib,1);
ExpRStar = sortedAllCalib(:,1);
FinalExpRStar = sortedAllCalib(:,2);
FinalIntCommands = sortedAllCalib(:,3);
indexNDF = sortedAllCalib(:,4);

%Ask for time of experiment, because circadian
prompt = {'Enter experiment time on a 24 hr clock'};
title = 'Experiment Time'; nLines = 1;
default = {'12'};
answer = inputdlg(prompt,title,nLines,default);
RecTime = str2double(answer{1});

%sizearray for seeing how the pupil changes
PupilChangeLow = zeros(length(ExpRStar),1);

smoothPupil = movmean(PupilTrace,SamplingRate/2,'omitnan');
[M, I] = min(smoothPupil((baseline*SamplingRate:(baseline+stimLength)*SamplingRate-SafetyFactor),:));
I = I+(baseline*SamplingRate);
%now calculate pupil area at specific times; for PupilChange5, it
%calculates the average of the Pupil area from 25s after stimulus onset to
%30s after stimulus onset
for ii = 1:length(ExpRStar)
    if short==1
        
        if (I(ii)+2.5*SamplingRate)>4188
            PupilChangeLow(ii) = mean(smoothPupil((I(ii)-2.5*SamplingRate):end,ii),'omitnan');
        else
            PupilChangeLow(ii) = mean(smoothPupil((I(ii)-2.5*SamplingRate):(I(ii)+2.5*SamplingRate),ii),'omitnan');
        end
    else
        PupilChangeLow(ii) = mean(smoothPupil((I(ii)-2.5*SamplingRate):(I(ii)+2.5*SamplingRate),ii),'omitnan');
    end
  
end

% Sigmoidal function for dose-response curve fitting
% Parameters: b(1)=minimum, b(2)=maximum, b(3)=EC50, b(4)=Hill slope
% Note: EC50 is the illuminance at halfway between min and max response
SigmoidFit = @(b,x)(b(1)+(b(2)-b(1))./(1+(b(3)./x).^b(4)));
%Initial guess for least squares regression
beta0 = [0,1,100,-1];
% Define bounds: minimum=0.2, maximum=1 (constrains to normalized range)
lb = [0.2,1,ExpRStar(1)/10,-10];  % Lower bounds
ub = [0.2,1,ExpRStar(end)*10,0];  % Upper bounds (Hill slope must be negative)




%Plot and Fit for PupilChangeLow
%subplot(5,1,1)
figure
points = semilogx(ExpRStar,PupilChangeLow, 'ko');
points.MarkerSize = 7;
points.MarkerFaceColor = [0 0 0];
ax = gca;
ax.Box = 'off';
ax.YLabel.String = 'Relative Pupil Area';
ax.YLabel.FontSize = 12;
ax.YLim = [0 1.2];
ax.YTick = [0.0 0.2 0.4 0.6 0.8 1.0];
ax.XLim = [ExpRStar(1)/10 ExpRStar(length(ExpRStar))*10];
ax.XTick = [10^0 10^1 10^2 10^3 10^4 10^5 10^6];
ax.XTickLabel = {'1', '10', '100', '1000', '10000', '100000', '1000000'};
ax.XLabel.String = 'Illuminance (R*)';
FitPupilChangeLow = lsqcurvefit(SigmoidFit,beta0,ExpRStar,PupilChangeLow,lb,ub);
EC50 = FitPupilChangeLow(3);
Hill = FitPupilChangeLow(4);
hold on
logrange = logspace(log10(ExpRStar(1)/10), log10(ExpRStar(length(ExpRStar))*10));
line = semilogx(logrange, SigmoidFit(FitPupilChangeLow,logrange));
line.Color = 'k';
line.LineWidth = 1;
maxref = refline(0,1);
maxref.Color = 'k';
maxref.LineWidth = 1;
maxref.LineStyle = '--';
minref = refline(0,0.2);
minref.Color = 'k';
minref.LineWidth = 1;
minref.LineStyle = '--';



%[savefileName, savepathName] = uiputfile('.mat','Save as');
savepathName = pathName;
names = strsplit(fileName,'_');
names(length(names)) = {'analyzed.mat'};
savefileName = strjoin(names,'_');
save([savepathName savefileName],'baseline', 'baselineAvg', 'darkIntensity','darkTime','ExpRStar', ...
    'FinalExpRStar','FinalIntCommands','indexNDF','NumRuns','param','PupilTrace', 'RecTime', 'RunsPerStim',...
    'SafetyFactor','SamplingRate','stimLength','TotalTime','EC50', 'Hill','PupilChangeLow','FitPupilChangeLow',...
    'smoothPupil','SigmoidFit','lb','ub')
disp('The End')
