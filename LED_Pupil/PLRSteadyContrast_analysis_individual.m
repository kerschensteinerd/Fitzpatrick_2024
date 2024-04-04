%This script analyzes PLR for a single experiment and generates a sigmoidal
%fit and curve. It then saves the values at each illuminance and the EC50
%from the curve fit. 
clear; clc; close all;

%%Load Parsed files
[fileName, pathName] = uigetfile('.mat', 'Load parsed file');
%Load first file to get parameters and size array
load([pathName fileName]);
disp(['Succesfully loaded ' fileName])


%Ask for time of experiment, because circadian
prompt = {'Enter experiment time on a 24 hr clock'};
title = 'Experiment Time'; nLines = 1;
default = {'12'};
answer = inputdlg(prompt,title,nLines,default);
RecTime = str2double(answer{1});

%sizearrays for seeing how the pupil changes

PupilChangeLow = zeros(NumRuns,1);

smoothPupil = movmean(PupilTrace,SamplingRate/2,'omitnan');
[M, I] = min(smoothPupil((baseline*SamplingRate:(baseline+stimLength)*SamplingRate),:));
I = I+(baseline*SamplingRate);


%now calculate pupil area at specific times; for PupilChange5, it
%calculates the average of the Pupil area from 5s after stimulus onset to
%10s after stimulus onset
for ii = 1:NumRuns
    
    PupilChangeLow(ii) = mean(smoothPupil((I(ii)-2.5*SamplingRate):(I(ii)+2.5*SamplingRate),ii),'omitnan');
    %figure
    %plot(PupilTrace(:,ii))
end


%Create matrices for each contrast

ContrastLow = zeros(length(ExpContrast)-1,length(frequency));

for ll = 2:NumRuns
    ContrastLow(mod(ll-2,length(ExpContrast)-1)+1,floor((ll+length(ExpContrast)-1-1.1)/(length(ExpContrast)-1))) = PupilChangeLow(ll);
end






%For each trace, generate a trace for analysis with mean = 0
cPupilTrace = PupilTrace(601:4188,:);
%ControlTrace = PupilTrace(1:600,:);
%ControlTrace = cPupilTrace;
FFTTrace = cPupilTrace;
for ii = 1:NumRuns
    FFTTrace(:,ii) = cPupilTrace(:,ii) - mean(cPupilTrace(:,ii), 'omitnan');
    
end
%Replace NaN with the mean
FFTTrace(isnan(FFTTrace)) = 0;


PupilFFT = fft(FFTTrace);
n = length(PupilFFT);
f = (0:n-1)*SamplingRate/n;

fshift = (-n/2:n/2-1)*SamplingRate/n;
PupilFFTShift = fftshift(PupilFFT,1);
PupilPower = abs(PupilFFTShift).^2/n;




%calculate index for each frequency
IFrequency = frequency;
for f = 1:length(frequency)
    [value, I] = min(abs(fshift-frequency(f)));
    IFrequency(f) = I;    
end


%Create matrices for each contrast
ContrastFFT = zeros(length(ExpContrast)-1,length(frequency));


for ll = 2:NumRuns
    ContrastFFT(mod(ll-2,length(ExpContrast)-1)+1,floor((ll+length(ExpContrast)-1-1.1)/(length(ExpContrast)-1))) = 2*PupilPower(IFrequency(floor((ll+length(ExpContrast)-1-1.1)/(length(ExpContrast)-1))), ll);
end


savepathName = pathName;
names = strsplit(fileName,'_');
names(length(names)) = {'analyzed.mat'};
savefileName = strjoin(names,'_');
save([savepathName savefileName],'baseline', 'baselineAvg', 'darkIntensity','darkTime','ExpRStar', ...
    'FinalRStarForContrast','FinalIntCommands','indexNDF','NumRuns','param','PupilTrace', 'smoothPupil','PupilChangeLow', 'PupilChange5', 'PupilChangeEnd','PupilChange5_early', 'PupilChange5_mid','RunsPerStim',...
    'SafetyFactor','SamplingRate','stimLength','TotalTime', 'Contrast', 'ContrastLow', 'ContrastEnd','ContrastMid', 'ContrastEarly', 'ContrastFFT', 'ExpContrast','frequency', 'RecTime')
disp('The End')