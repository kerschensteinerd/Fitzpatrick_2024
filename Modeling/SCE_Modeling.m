%% Stiles-Crawford Effect (SCE) Modeling
% This script models the Stiles-Crawford Effect and calculates retinal
% illuminance as a function of pupil size, both with and without SCE.
% It also generates optical transfer functions for different photoreceptor types.
%
% The Stiles-Crawford Effect describes the reduced efficiency of light
% entering the eye at the pupil periphery compared to the center.
%
% Requirements:
%   - MATLAB Image Processing Toolbox
%   - isetbio toolbox (https://github.com/isetbio/isetbio)
%   - mouseCore function (custom - contact authors)
%
% Author: Fitzpatrick et al., 2024

%% initialize
clear; clc; close all;

% ========== CONFIGURATION - UPDATE THIS PATH ==========
% Path to isetbio toolbox installation
% Download from: https://github.com/isetbio/isetbio
isetbio_path = 'X:\\isetbio-master';  % UPDATE THIS PATH
if ~exist(isetbio_path, 'dir')
    warning(['isetbio path not found: ' isetbio_path]);
    warning('Please update isetbio_path variable to your isetbio installation directory');
    warning('Download from: https://github.com/isetbio/isetbio');
    error('Cannot proceed without isetbio');
end
addpath(genpath(isetbio_path))
% ========== END CONFIGURATION ==========

%import normalized photoreceptor spectra (area = 1)
%corrected for pre-receptoral filtering
load SCone
load ipRGC
load Rod
load MCone


ieInit;

%% Parameters
% Wavelength range for spectral analysis
StartWaveLen = 300;  % nm
EndWaveLen = 780;    % nm
StepWaveLen = 5;     % nm

% Simulation parameters
nCSF = 500;     % Number of points for contrast sensitivity function
nPupil = 100;   % Number of pupil sizes to test
nIllum = 100;   % Number of illuminance levels
nFreq = 500;    % Number of spatial frequencies
nDist = 500;    % Number of distances

% Mouse eye optical parameters (from Remtulla & Hallett)
fl = 2.347/1000; %focal length in mm, back half of mouse eye
tr = 0.7;% transmittance of the ocular media

% Stiles-Crawford Effect parameters
% Based on FWHM estimated from Wei Li's work in ground squirrel
FWHM_d = 12; %in degrees, full width of half maximum
FWHM_m = fl*tand(FWHM_d);  % Convert to meters
SCE_sigma = FWHM_m/(2*sqrt(2*log(2)));  % Gaussian sigma from FWHM

%% calculated params

% Illumination and contrast parameters
IllumArea = 1.02;  % Illumination area in mm²
IllumR = sqrt(IllumArea/pi)/1000;  % Illumination radius in meters
RelF0 = 0.73;  % Relative area for contrast condition
ContR = sqrt(RelF0/pi)/1000;  % Contrast radius in meters

% Spatial frequency range (log-spaced)
freqSf = logspace(-1.5,1,nFreq);  % cycles per degree

% Distance range for depth of field calculations
dof_dist = linspace(0.1/100,50/100,nDist);  % 0.1 to 50 cm in meters

% Dioptric power calculation
mD0 = 1/fl*1.3341; %diopteric power = (1/focal length in m) * refractive index of vitreous

% Pupil radius range
mPupilRadius = linspace(0.3,0.85,nPupil)/1000; %in meters (0.3-0.85 mm)

% Illuminance range (log-spaced from 1 to 100,000 R*)
Illum = logspace(0,5,nIllum);

% 2D Gaussian function for Stiles-Crawford Effect
gaus2d = @(x,y,xmu,ymu,sig)exp(-((((x-xmu).^2)/(2*sig.^2)) + (((y-ymu).^2)/(2*sig.^2))));

%% Illuminance calculation, with and without Stiles Crawford Effect

Illum_Mat = zeros(nPupil,nIllum);
SCE_Mat = zeros(nPupil,nIllum);
for i = 1:nPupil
    Illum_Mat(i,:) = (Illum*tr*((mPupilRadius(i))^2)/(2*fl^2));
    %modified from Sliney et al. 2002 and Lyubarsky 2004, using illuminance on both sides
    %int
    
    Pup_x = linspace(-mPupilRadius(i),mPupilRadius(i),100);
    Pup_y = linspace(-mPupilRadius(i),mPupilRadius(i),100);
    [SCE_x,SCE_y] = meshgrid(Pup_x,Pup_y);
    SCE_gauss = gaus2d(SCE_x,SCE_y,0,0,SCE_sigma);
    
    SCE_Mat(i,:) = mean(SCE_gauss,'all') .* Illum_Mat(i,:);
end

[plotIllum, plotPupil1] = meshgrid(Illum,mPupilRadius);

rel_Illum_Mat = Illum_Mat./Illum_Mat(end,:);
rel_SCE_Mat = SCE_Mat./SCE_Mat(end,:);



%start here, reworking these to be flipped the other way
figure
subplot(2,1,1)
Illum_fig = pcolor(plotPupil1'.*1000,log10(plotIllum'),log10(Illum_Mat'));
Illum_fig.EdgeColor = 'none';
ax = gca;
ax.YLabel.String = 'Environmental Illuminance (log10 R*)';
ax.XLabel.String = 'Pupil Radius (mm)';
caxis([-1,4])
c = colorbar();
c.Label.String = 'Retinal Illuminance (log10 R*)';

subplot(2,1,2)
SCE_fig = pcolor(plotPupil1'.*1000,log10(plotIllum'),log10(SCE_Mat'));
SCE_fig.EdgeColor = 'none';
ax = gca;
ax.YLabel.String = 'Environmental Illuminance (log10 R*)';
ax.XLabel.String = 'Pupil Radius (mm)';
caxis([-1,4])
c = colorbar();
c.Label.String = 'Retinal Illuminance (log10 R*)';


%% mouse OTF

wave = StartWaveLen:StepWaveLen:EndWaveLen;

%mouse otf as a function of pupil radius for different
%photoreceptor-weighted wavelengths
S_Mat = zeros(nPupil,nFreq);
M_Mat = zeros(nPupil,nFreq);
R_Mat = zeros(nPupil,nFreq);
ip_Mat = zeros(nPupil,nFreq);
for i = 1:length(mPupilRadius)

    [otf, ~] = mouseCore(wave,freqSf, mPupilRadius(i), mD0);
    
    s_otf  = sum(otf.*SCone(:,2),1);
    S_Mat(i,:) = s_otf;
    
    
    m_otf  = sum(otf.*MCone(:,2),1);
    M_Mat(i,:) = m_otf;
  
    r_otf  = sum(otf.*Rod(:,2),1);
    R_Mat(i,:) = r_otf;
   
    ip_otf  = sum(otf.*ipRGC(:,2),1);
    ip_Mat(i,:) = ip_otf;
    
end

[plotSF, plotPupil] = meshgrid(freqSf,mPupilRadius);

figure
subplot(4,1,1)
S_fig = pcolor(plotPupil'.*1000,log10(plotSF'),abs(S_Mat'));
S_fig.EdgeColor = 'none';
ax = gca;
ax.YLabel.String = 'Spatial Frequency (log10 cpd)';
ax.XLabel.String = 'Pupil Radius (mm)';
caxis([0,1])
c = colorbar();
c.Label.String = 'Percent Maximum';

subplot(4,1,2)
M_fig = pcolor(plotPupil'.*1000,log10(plotSF'),abs(M_Mat'));
M_fig.EdgeColor = 'none';
ax = gca;
ax.YLabel.String = 'Spatial Frequency (log10 cpd)';
ax.XLabel.String = 'Pupil Radius (mm)';
caxis([0,1])
c = colorbar();
c.Label.String = 'Percent Maximum';

subplot(4,1,3)
R_fig = pcolor(plotPupil'.*1000,log10(plotSF'),abs(R_Mat'));
R_fig.EdgeColor = 'none';
ax = gca;
ax.YLabel.String = 'Spatial Frequency (log10 cpd)';
ax.XLabel.String = 'Pupil Radius (mm)';
caxis([0,1])
c = colorbar();
c.Label.String = 'Percent Maximum';


subplot(4,1,4)
ip_fig = pcolor(plotPupil'.*1000,log10(plotSF'),abs(ip_Mat'));
ip_fig.EdgeColor = 'none';
ax = gca;
ax.YLabel.String = 'Spatial Frequency (log10 cpd)';
ax.XLabel.String = 'Pupil Radius (mm)';
caxis([0,1])
c = colorbar();
c.Label.String = 'Percent Maximum';

