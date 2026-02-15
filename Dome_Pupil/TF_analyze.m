%% Temporal Frequency (TF) Analysis - Dome Setup
% This script analyzes pupillary responses to temporally modulated stimuli
% at different spatial and temporal frequencies. It performs FFT analysis
% to extract power at specific stimulus frequencies.
%
% Input files (selected via GUI):
%   - Parameter file: Contains stimIn/stimOut with spatial/temporal frequency parameters
%   - Left eye trace: Pupil area measurements over time
%   - Right eye trace: Pupil area measurements over time
%
% Outputs:
%   - FFT power at stimulus frequencies
%   - Pupil area organized by spatial frequency, temporal frequency, and repeat
%   - Saved .mat file with all analysis results
%
% Author: Fitzpatrick et al., 2024

clear; clc; close all;

%% Define variables and load data

% Experimental parameters - adjust these for your setup
stim_fps = 60;  % Stimulus presentation frame rate (Hz)
movie_fps = 15; % Pupil tracking camera frame rate (Hz)
mm_per_pix = 0.00345; % Conversion: millimeters per pixel for pupil measurements

% Analysis options
fixed = 0; % 0: use dynamic window around minimum, 1: use fixed 5s window
Global = 1; % Global analysis flag

%locomotion filtering thresholds
loco_filter = 1; % 1: enable locomotion filtering, 0: disable
spd_thresh = 0.5; % Speed threshold in cm/s to classify as locomotion
t_thresh = 2; % Minimum duration of running bout (seconds)
extra = 30; % Additional time after running to exclude (seconds) for pupil recovery

%parameters
[fileName, pathName]=uigetfile ('*.mat','Select parameter data');
load([pathName fileName]);

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

% Extract speed/distance of mouse from wheel encoder
% Wheel parameters: 15 cm diameter, 4028 pulses per rotation
mSpeed = zeros(length(Data),1);
mSpeed(2:end) = abs(diff(Data(:,8)));
mSpeed(mSpeed>150 | mSpeed<-150)= 0; % Remove points where encoder cycles (bit limit)
% Convert to cm/s: (pulses/frame) * (circumference/pulses_per_rotation) * fps
mSpeed = mSpeed*(15*pi*stim_fps/4028); % now in cm/s
mSpeed = downsample(mSpeed,stim_fps/movie_fps); % Match pupil tracking rate
smoothSpeed = movmean(mSpeed,movie_fps/2,'omitnan'); % Smooth with 0.5s window

% Apply locomotion filtering if enabled
if loco_filter
    BW = smoothSpeed>spd_thresh;  % Identify frames above threshold
    locofiltID = bwareaopen(BW,t_thresh*movie_fps);  % Keep only sustained bouts
    % Extend blanked region to account for post-locomotion pupil recovery
    [labeledImage, numBlobs] = bwlabel(locofiltID);
    measurements = regionprops(labeledImage, 'BoundingBox');
    for k = 1 : numBlobs
        % Get the bounding box for this locomotion bout
        thisBB = measurements(k).BoundingBox;
        % Get rows of the bounding box
        row1 = ceil(thisBB(2));
        row2 = row1 + thisBB(4)+extra*movie_fps;  % Add recovery time
        if row2>length(LArea)
            row2 = length(LArea);
        end
        % Mark these frames as invalid
        locofiltID(row1:row2) = 1;
    end
    LArea(locofiltID) = NaN;  % Set locomotion frames to NaN
    RArea(locofiltID) = NaN;
end



%% Process data
time = Data(:,1);

cpd = stimIn.cpd;
nCpd = length(cpd);
cpdOrder = stimOut.stimcpd;

freq = stimIn.freq;
nFreq = length(freq);
freqOrder = stimOut.stimfreq;

nRepeats =stimIn.nRepeats;

nStim = nRepeats*nFreq*nCpd;

%script cycles in this order --> cpd, then freq, then repeats 
cpdIdx = repmat(stimOut.stimcpd,1,length(stimIn.freq)*stimIn.nRepeats);
freqIdx = double(repmat(repelem(stimOut.stimfreq,length(stimIn.cpd)),1,stimIn.nRepeats));
nRepeatIdx = double(repelem(1:stimIn.nRepeats,1,length(stimIn.freq)*length(stimIn.cpd)));



runTime = double(stimIn.interDur+stimIn.Pre+stimIn.Dur+stimIn.Post)*stim_fps;
PupTime = double(runTime/stim_fps*movie_fps);

Speed_runs = reshape(smoothSpeed,PupTime,length(stimIn.cpd)*length(stimIn.freq)*stimIn.nRepeats);


%% Calculate Pupil Traces

LArea_runs = reshape(LArea,PupTime,length(stimIn.cpd)*length(stimIn.freq)*stimIn.nRepeats);
LArea_pre = mean(LArea_runs(stimIn.interDur*movie_fps:(stimIn.interDur+stimIn.Pre)*movie_fps,:),'omitnan');
LArea_pre = reshape(LArea_pre,length(stimIn.cpd),length(stimIn.freq),stimIn.nRepeats);

RArea_runs = reshape(RArea,PupTime,length(stimIn.cpd)*length(stimIn.freq)*stimIn.nRepeats);
RArea_pre = mean(RArea_runs(stimIn.interDur*movie_fps:(stimIn.interDur+stimIn.Pre)*movie_fps,:),'omitnan');
RArea_pre = reshape(RArea_pre,length(stimIn.cpd),length(stimIn.freq),stimIn.nRepeats);

smoothLPupil = movmean(LArea_runs,movie_fps/2,'omitnan');
[L_M, L_I] = min(smoothLPupil(((stimIn.interDur+stimIn.Pre+stimIn.Dur/2)*movie_fps:(stimIn.interDur+stimIn.Pre+stimIn.Dur)*movie_fps),:));
%L_nanidx = sum(isnan(smoothLPupil(((stimIn.interDur+stimIn.Pre+stimIn.Dur/2)*movie_fps:(stimIn.interDur+stimIn.Pre+stimIn.Dur)*movie_fps),:)));
L_I = L_I+double((stimIn.interDur+stimIn.Pre+stimIn.Dur/2)*movie_fps);
%L_I(L_nanidx>round(stimIn.Dur/4)) = NaN;

smoothRPupil = movmean(RArea_runs,movie_fps/2,'omitnan');
[R_M, R_I] = min(smoothRPupil(((stimIn.interDur+stimIn.Pre+stimIn.Dur/2)*movie_fps:(stimIn.interDur+stimIn.Pre+stimIn.Dur)*movie_fps),:));
%R_nanidx = sum(isnan(smoothRPupil(((stimIn.interDur+stimIn.Pre+stimIn.Dur/2)*movie_fps:(stimIn.interDur+stimIn.Pre+stimIn.Dur)*movie_fps),:)));
R_I = R_I+double((stimIn.interDur+stimIn.Pre+stimIn.Dur/2)*movie_fps);
%R_I(R_nanidx>round(stimIn.Dur/4)) = NaN;



%% Calculate power

FFTLPupil_trace = smoothLPupil-mean(smoothLPupil,'omitnan');
FFTLPupil_trace(isnan(FFTLPupil_trace)) = 0;
FFTRPupil_trace = smoothRPupil-mean(smoothRPupil,'omitnan');
FFTRPupil_trace(isnan(FFTRPupil_trace)) = 0;

LPupFFT = fft(FFTLPupil_trace);
RPupFFT = fft(FFTRPupil_trace);


n = double(PupTime);
fshift = (-n/2:n/2-1)*movie_fps/n;

LPupFFTShift = fftshift(LPupFFT,1);
RPupFFTShift = fftshift(RPupFFT,1);
LPupFFTPower = abs(LPupFFTShift).^2/n;
RPupFFTPower = abs(RPupFFTShift).^2/n;


%calculate index for each frequency
LPupPower = zeros(1,length(stimIn.cpd)*length(stimIn.freq)*stimIn.nRepeats);
RPupPower = zeros(1,length(stimIn.cpd)*length(stimIn.freq)*stimIn.nRepeats);
for f = 1:length(freqIdx)
    if freqIdx(f)>max(fshift)
        LPupPower(f) = NaN;
        RPupPower(f) = NaN;
    else
        [value, I] = min(abs(fshift-freqIdx(f)));
        LPupPower(f) = 2*max(LPupFFTPower(I-1:I+1,f));
        RPupPower(f) = 2*max(RPupFFTPower(I-1:I+1,f));
    end
end

LPupPower = reshape(LPupPower,length(stimIn.cpd),length(stimIn.freq),stimIn.nRepeats);
RPupPower = reshape(RPupPower,length(stimIn.cpd),length(stimIn.freq),stimIn.nRepeats);


%% distribute other variables

PupilChangeLeft = zeros(1,length(stimIn.cpd)*length(stimIn.freq)*stimIn.nRepeats);
PupilChangeRight = zeros(1,length(stimIn.cpd)*length(stimIn.freq)*stimIn.nRepeats);

RPup = zeros(length(stimIn.cpd),length(stimIn.freq),stimIn.nRepeats);
LPup = zeros(length(stimIn.cpd),length(stimIn.freq),stimIn.nRepeats);

RPupTrace = zeros(length(stimIn.cpd),length(stimIn.freq),stimIn.nRepeats,length(smoothRPupil));
LPupTrace = zeros(length(stimIn.cpd),length(stimIn.freq),stimIn.nRepeats,length(smoothLPupil));
SpeedTrace = zeros(length(stimIn.cpd),length(stimIn.freq),stimIn.nRepeats,length(Speed_runs));

for i = 1:length(stimIn.cpd)*length(stimIn.freq)*stimIn.nRepeats
    
    if fixed==1
        PupilChangeLeft(i) = mean(smoothLPupil((stimIn.interDur+stimIn.Pre+stimIn.Dur/2)*movie_fps:(stimIn.interDur+stimIn.Pre+stimIn.Dur)*movie_fps,i),'omitnan');
        PupilChangeRight(i) = mean(smoothRPupil((stimIn.interDur+stimIn.Pre+stimIn.Dur/2)*movie_fps:(stimIn.interDur+stimIn.Pre+stimIn.Dur)*movie_fps,i),'omitnan');
    else
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
    [c,f,n] = ind2sub([length(stimIn.cpd),length(stimIn.freq),stimIn.nRepeats],i);
    RPup(c,f,n) = PupilChangeRight(i);
    RPupTrace(c,f,n,:) = smoothRPupil(:,i);
    LPup(c,f,n) = PupilChangeLeft(i);
    LPupTrace(c,f,n,:) = smoothLPupil(:,i);
    SpeedTrace(c,f,n,:) = Speed_runs(:,i);
end


%% sort and save
[cpd, I_cpd] = sort(stimOut.stimcpd);
[freq, I_freq] = sort(stimOut.stimfreq);



RPup = RPup(I_cpd,I_freq,:);
LPup = LPup(I_cpd,I_freq,:);
RArea_pre = RArea_pre(I_cpd,I_freq,:);
LArea_pre = LArea_pre(I_cpd,I_freq,:);

RPupPower = RPupPower(I_cpd,I_freq,:);
LPupPower = LPupPower(I_cpd,I_freq,:);
RPupTrace = RPupTrace(I_cpd,I_freq,:,:);
LPupTrace = LPupTrace(I_cpd,I_freq,:,:);

%SpeedTrace = SpeedTrace(I_cpd,I_freq,:,:);

save_pathName = pathName;
names = strsplit(fileName,'_');
names(length(names)) = {'analyzed.mat'};
save_fileName = strjoin(names,'_');
save([save_pathName save_fileName]);

