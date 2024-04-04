%% initialize
clear; clc; close all;

addpath(genpath('X:\\isetbio-master'))

%import normalized photoreceptor spectra (area = 1)
%corrected for pre-receptoral filtering
load SCone
load ipRGC
load Rod
load MCone


ieInit;

%% Parameters
StartWaveLen = 300;
EndWaveLen = 780;
StepWaveLen = 5;



freqSf = logspace(-1.5,1,500);

mD0 = 1/0.001756; %diopteric power of mouse eye = 1/focal length
mPupilRadius = linspace(0.3,0.9,10)/1000; %in m


IllumArea = 1.02;
IllumR = sqrt(IllumArea/pi)/1000;
RelF0 = 0.73;
ContR = sqrt(RelF0/pi)/1000;

%video settings

video = 'X:\\NaturalMovies\\RawData\\CT9000\\Nature\\All\\20190616_180258_3.mp4';
frame = 1;

fullHeight = 1080;
fullWidth = 1920;

%visual angle to clip to
ClipAngH = 50;
ClipAngW = 50;

CamAngH = 55; % angle
CamAngW = 100; % angle

vHeight = round(fullHeight*ClipAngH/CamAngH);
vWidth = round(fullWidth*ClipAngW/CamAngW);

centeredHeight = ((fullHeight-vHeight)/2)+1:fullHeight-((fullHeight-vHeight)/2);
centeredWidth = ((fullWidth-vWidth)/2)+1:fullWidth-((fullWidth-vWidth)/2);

%% mouse OTF

wave = StartWaveLen:StepWaveLen:EndWaveLen;

%mouse otf as a function of pupil radius for different
%photoreceptor-weighted wavelengths
figure;
for i = 1:length(mPupilRadius)

    [otf, ~] = mouseCore(wave,freqSf, mPupilRadius(i), mD0);
    
    subplot(4,2,1)
    hold on
    s_otf  = sum(otf.*SCone(:,2),1);
    plot(log10(freqSf),s_otf,'Color',repmat((i-1)/length(mPupilRadius),[1,3]))
    
    subplot(4,2,3)
    hold on
    m_otf  = sum(otf.*MCone(:,2),1);
    plot(log10(freqSf),m_otf,'Color',repmat((i-1)/length(mPupilRadius),[1,3]))
    
    subplot(4,2,5)
    hold on
    r_otf  = sum(otf.*Rod(:,2),1);
    plot(log10(freqSf),r_otf,'Color',repmat((i-1)/length(mPupilRadius),[1,3]))
    
    subplot(4,2,7)
    hold on
    ip_otf  = sum(otf.*ipRGC(:,2),1);
    plot(log10(freqSf),ip_otf,'Color',repmat((i-1)/length(mPupilRadius),[1,3]))
    
end

%now plotted specifically for illumR and contR
%figure;
[i_otf, ~] = mouseCore(wave,freqSf, IllumR, mD0);
[c_otf, ~] = mouseCore(wave,freqSf, ContR, mD0);

subplot(4,2,2)
s_i_otf  = sum(i_otf.*SCone(:,2),1);
plot(log10(freqSf),s_i_otf,'g')
hold on
s_c_otf  = sum(c_otf.*SCone(:,2),1);
plot(log10(freqSf),s_c_otf,'k')

subplot(4,2,4)
m_i_otf  = sum(i_otf.*MCone(:,2),1);
plot(log10(freqSf),m_i_otf,'g')
hold on
m_c_otf  = sum(c_otf.*MCone(:,2),1);
plot(log10(freqSf),m_c_otf,'k')

subplot(4,2,6)
r_i_otf  = sum(i_otf.*Rod(:,2),1);
plot(log10(freqSf),r_i_otf,'g')
hold on
r_c_otf  = sum(c_otf.*Rod(:,2),1);
plot(log10(freqSf),r_c_otf,'k')

subplot(4,2,8)
ip_i_otf  = sum(i_otf.*ipRGC(:,2),1);
plot(log10(freqSf),ip_i_otf,'g')
hold on
ip_c_otf  = sum(c_otf.*ipRGC(:,2),1);
plot(log10(freqSf),ip_c_otf,'k')


figure
for i = 1:length(mPupilRadius)

    [otf, ~] = mouseCore(wave,freqSf, mPupilRadius(i), mD0);
    
    subplot(4,1,1)
    hold on
    s_otf  = sum(otf.*SCone(:,2),1);
    plot(log10(freqSf),s_otf,'Color',repmat((i-1)/length(mPupilRadius),[1,3]))
    
    subplot(4,1,2)
    hold on
    m_otf  = sum(otf.*MCone(:,2),1);
    plot(log10(freqSf),m_otf,'Color',repmat((i-1)/length(mPupilRadius),[1,3]))
    
    subplot(4,1,3)
    hold on
    r_otf  = sum(otf.*Rod(:,2),1);
    plot(log10(freqSf),r_otf,'Color',repmat((i-1)/length(mPupilRadius),[1,3]))
    
    subplot(4,1,4)
    hold on
    ip_otf  = sum(otf.*ipRGC(:,2),1);
    plot(log10(freqSf),ip_otf,'Color',repmat((i-1)/length(mPupilRadius),[1,3]))
    
end

%% now image analysis

v=VideoReader(video);
f = read(v,frame);
mono_f =mean(f(centeredHeight,centeredWidth,2:3),3)/255; %make monochromatic of G and B

figure
imshow(imrotate(mono_f,180),[0 1]);

%linearize image
lin_f = zeros(size(mono_f));
for i = 1:numel(mono_f)
    if mono_f(i)<=0.04045
        lin_f(i) = mono_f(i)/12.92;
    else
        lin_f(i) = ((mono_f(i) + 0.055)/(1.055))^2.4;
    end
end

%fft of image
fft_f = fft2(lin_f);
fft_fshift =fftshift(fft_f);

SFx = vHeight/ClipAngH; %pixels per degree = spatial sampling rate
x = (0:vHeight-1) * (SFx/vHeight);
xl = x(x<SFx/2);
nx = length(xl);

SFy = vWidth/ClipAngW; %pixels per degree = spatial sampling rate
y = (0:vWidth-1) * (SFy/vWidth);
yl = y(y<SFy/2);
ny = length(yl);



%imagesc(log10(abs(fftshift(fft_f)).^2/(vWidth*vHeight)))

%can't use 2d otf, because it forces bad spatial resolution to be fast
%create it myself from mouse core
%use x because longer dimension in sf space
[i_otfx, ~] = mouseCore(wave,xl, IllumR, mD0);
[c_otfx, ~] = mouseCore(wave,xl, ContR, mD0);
[d_otfx, ~] = mouseCore(wave,xl, mPupilRadius(end), mD0);

%SCone
%illum 
subplot(4,3,1)

s_i_otfx  = sum(i_otfx.*SCone(:,2),1);
[s_i_mono] = mSee(fft_fshift,s_i_otfx,xl,yl);
imshow(imrotate(s_i_mono,180),[0 1]);

%contrast 
subplot(4,3,2)

s_c_otfx  = sum(c_otfx.*SCone(:,2),1);
[s_c_mono] = mSee(fft_fshift,s_c_otfx,xl,yl);
imshow(imrotate(s_c_mono,180),[0 1]);

%delta MTF
subplot(4,3,3)
delta_s_otfx = s_c_otfx-s_i_otfx;
plot(log10(xl),delta_s_otfx)
ax = gca;
ax.XLim = [-1.5 1];
ax.YLim = [-0.15 .25];
hold on
yline(0,'--k')

%MCone
%illum 
subplot(4,3,4)

m_i_otfx  = sum(i_otfx.*MCone(:,2),1);
[m_i_mono] = mSee(fft_fshift,m_i_otfx,xl,yl);
imshow(imrotate(m_i_mono,180),[0 1]);

%contrast 
subplot(4,3,5)

m_c_otfx  = sum(c_otfx.*MCone(:,2),1);
[m_c_mono] = mSee(fft_fshift,m_c_otfx,xl,yl);
imshow(imrotate(m_c_mono,180),[0 1]);

m_d_otfx = sum(d_otfx.*MCone(:,2),1);
[m_d_mono] = mSee(fft_fshift,m_d_otfx,xl,yl);
imshow(imrotate(m_d_mono,180),[0 1]);

%delta MTF
subplot(4,3,6)
delta_m_otfx = m_c_otfx-m_i_otfx;
plot(log10(xl),delta_m_otfx)
ax = gca;
ax.XLim = [-1.5 1];
ax.YLim = [-0.15 .25];
hold on
yline(0,'--k')

%Rod
%illum 
subplot(4,3,7)

r_i_otfx  = sum(i_otfx.*Rod(:,2),1);
[r_i_mono] = mSee(fft_fshift,r_i_otfx,xl,yl);
imshow(imrotate(r_i_mono,180),[0 1]);

%contrast 
subplot(4,3,8)

r_c_otfx  = sum(c_otfx.*Rod(:,2),1);
[r_c_mono] = mSee(fft_fshift,r_c_otfx,xl,yl);
imshow(imrotate(r_c_mono,180),[0 1]);

%delta MTF
subplot(4,3,9)
delta_r_otfx = r_c_otfx-r_i_otfx;
plot(log10(xl),delta_r_otfx)
ax = gca;
ax.XLim = [-1.5 1];
ax.YLim = [-0.15 .25];
hold on
yline(0,'--k')

%iprgc
%illum 
subplot(4,3,10)

ip_i_otfx  = sum(i_otfx.*ipRGC(:,2),1);
[ip_i_mono] = mSee(fft_fshift,ip_i_otfx,xl,yl);
imshow(imrotate(ip_i_mono,180),[0 1]);

%contrast 
subplot(4,3,11)

ip_c_otfx  = sum(c_otfx.*ipRGC(:,2),1);
[ip_c_mono] = mSee(fft_fshift,ip_c_otfx,xl,yl);
imshow(imrotate(ip_c_mono,180),[0 1]);

%delta MTF
subplot(4,3,12)
delta_ip_otfx = ip_c_otfx-ip_i_otfx;
plot(log10(xl),delta_ip_otfx)
ax = gca;
ax.XLim = [-1.5 1];
ax.YLim = [-0.15 .25];
hold on
yline(0,'--k')



