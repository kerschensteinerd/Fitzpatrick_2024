# Pupillary Contrast Response Analysis - Fitzpatrick et al. 2024

This repository contains MATLAB analysis and modeling code supporting the manuscript:

**"A pupillary contrast response in mice and humans - Neural mechanisms and visual functions"**  
*Fitzpatrick et al., 2024*

## Overview

This codebase provides tools for analyzing pupillary light reflex (PLR) and pupillary contrast responses in mice and humans, along with computational models of mouse vision through different photoreceptor types.

### Key Features

- **Pupil tracking analysis** with area measurements and constriction quantification
- **Sigmoidal curve fitting** for dose-response relationships (Hill equation, EC50 calculation)
- **Temporal frequency analysis** using FFT-based power spectral analysis
- **Locomotion filtering** to remove behavioral artifacts during mouse running
- **Optical transfer function (OTF) modeling** for mouse photoreceptors (S-cone, M-cone, Rod, ipRGC)
- **Stiles-Crawford Effect simulation** for retinal illuminance modeling

## Repository Structure

```
.
├── Dome_Pupil/                  # Analysis for projection-based experiments
│   ├── PLR_analyze.m           # Pupillary light reflex analysis with sigmoidal fitting
│   └── TF_analyze.m            # Temporal frequency response analysis
├── LED_Pupil/                   # Analysis for LED-based experiments
│   ├── PLR_analysis_individual.m              # Individual PLR analysis
│   └── PLRSteadyContrast_analysis_individual.m # Contrast response analysis
├── Modeling/                    # Computational vision models
│   ├── MF_Modeling.m           # Multi-photoreceptor OTF modeling
│   ├── SCE_Modeling.m          # Stiles-Crawford Effect modeling
│   ├── mSee.m                  # OTF filter application utility
│   └── Data/                   # Photoreceptor spectral sensitivity data
│       ├── SCone.mat
│       ├── MCone.mat
│       ├── Rod.mat
│       ├── ipRGC.mat
│       └── NatScene.png
└── Fitzpatrick et al. 2024 - [...].pdf  # Published manuscript
```

## Requirements

### Software Requirements

- **MATLAB** (tested with R2019b or later)
  - Image Processing Toolbox
  - Optimization Toolbox (for `lsqcurvefit`)
  - Signal Processing Toolbox (for FFT operations)

### External Dependencies

#### For Modeling Scripts Only

The modeling scripts (`MF_Modeling.m` and `SCE_Modeling.m`) require:

- **[isetbio](https://github.com/isetbio/isetbio)** - Image Systems Engineering Toolbox for Biology
  - Installation: Clone or download from https://github.com/isetbio/isetbio
  - Configure the path in the scripts (see Configuration section below)

- **mouseCore function** - Custom optical modeling function (not included)
  - Contact the authors if you need this function for OTF calculations

## Installation

1. **Clone this repository:**
   ```bash
   git clone https://github.com/kerschensteinerd/Fitzpatrick2024_pupillary-contrast-response.git
   cd Fitzpatrick2024_pupillary-contrast-response
   ```

2. **Set up MATLAB environment:**
   - Ensure you have MATLAB with required toolboxes installed
   - Add this repository to your MATLAB path

3. **Install isetbio (for modeling only):**
   ```bash
   git clone https://github.com/isetbio/isetbio.git
   ```
   - Update the path in modeling scripts (see Configuration below)

## Configuration

### Configuring Paths for Modeling Scripts

Before running `MF_Modeling.m` or `SCE_Modeling.m`, update the following paths at the top of each file:

```matlab
% Update this path to your isetbio installation
addpath(genpath('/path/to/your/isetbio-master'))

% For MF_Modeling.m only - update to your video file location
video = '/path/to/your/natural/movies/video.mp4';
```

### Experimental Parameters

Key parameters that may need adjustment for your experiments:

**Dome_Pupil/PLR_analyze.m:**
- `stim_fps = 60;` - Stimulus frame rate
- `movie_fps = 15;` - Pupil tracking frame rate
- `mm_per_pix = 0.00345;` - Conversion factor for pupil measurements
- `cd_to_R = 13.8893;` - Empirical conversion from cd/m² to R* (photoisomerizations)

**LED_Pupil scripts:**
- `loco_filter = 1;` - Enable/disable locomotion filtering
- `spd_thresh = 0.5;` - Speed threshold (cm/s) for locomotion filtering
- `t_thresh = 2;` - Minimum running bout duration (s)

## Usage

### 1. Pupillary Light Reflex Analysis (Dome Setup)

**Purpose:** Analyze pupil area changes in response to different illuminance levels and fit sigmoidal dose-response curves.

```matlab
% Open MATLAB and navigate to the Dome_Pupil directory
cd Dome_Pupil

% Run the analysis script
PLR_analyze

% Follow GUI prompts to:
% 1. Select parameter data file (.mat)
% 2. Select left eye trace file (.mat)
% 3. Select right eye trace file (.mat)
% 4. Choose save location for results
```

**Outputs:**
- Sigmoidal curve plots for left and right eyes
- EC50 values (half-maximal effective concentration)
- Hill slope coefficients
- Saved `.mat` file with all analysis results

### 2. Temporal Frequency Analysis (Dome Setup)

**Purpose:** Analyze pupil responses to temporally modulated stimuli at different spatial frequencies.

```matlab
cd Dome_Pupil
TF_analyze

% Follow GUI prompts to select data files
```

**Outputs:**
- FFT power spectra for each stimulus condition
- Pupil area traces organized by spatial frequency and temporal frequency
- FFT power values at stimulus frequencies

### 3. Individual PLR Analysis (LED Setup)

**Purpose:** Analyze individual subject pupillary responses to LED stimuli.

```matlab
cd LED_Pupil
PLR_analysis_individual

% Provide experiment time when prompted (24-hour format)
```

**Outputs:**
- Sigmoidal dose-response curves
- EC50 and Hill coefficients
- Smoothed pupil traces

### 4. Steady Contrast Analysis (LED Setup)

**Purpose:** Analyze pupil responses to contrast stimuli with FFT analysis.

```matlab
cd LED_Pupil
PLRSteadyContrast_analysis_individual

% Provide experiment time when prompted
```

**Outputs:**
- Contrast response matrices
- FFT power analysis for different frequencies and contrasts

### 5. Optical Transfer Function Modeling

**Purpose:** Model how pupil size affects spatial frequency transmission through mouse photoreceptors.

```matlab
cd Modeling

% For multi-photoreceptor OTF modeling and natural scene filtering:
MF_Modeling

% For Stiles-Crawford Effect modeling:
SCE_Modeling
```

**Note:** These require isetbio and the mouseCore function (see Requirements).

**Outputs:**
- OTF plots as functions of pupil radius and spatial frequency
- Filtered natural scene images for different photoreceptor types
- Retinal illuminance calculations with and without Stiles-Crawford Effect

## Data Format

### Input Files

All analysis scripts expect MATLAB `.mat` files containing:

**Parameter files:**
- `stimIn` - Structure with stimulus parameters (duration, illuminance, frequencies, etc.)
- `stimOut` - Array with stimulus presentation order
- `Data` - Trial-by-trial data matrix
- `corrData` - Locomotion/wheel speed data (for Dome_Pupil)

**Trace files:**
- `area` - Vector of pupil area measurements over time (in pixels²)

### Output Files

Analysis scripts save comprehensive `.mat` files including:
- Original input parameters
- Processed pupil traces
- Fitted curve parameters (EC50, Hill slope)
- Statistical measures
- FFT results (for frequency analysis)

## Key Parameters and Constants

### Conversion Factors

- **mm_per_pix = 0.00345** - Millimeters per pixel for pupil area measurements (camera-dependent)
- **cd_to_R = 13.8893** - Conversion from cd/m² to R* (rod photoisomerizations per photoreceptor per second)

### Locomotion Filtering

- **spd_thresh = 0.5 cm/s** - Minimum speed to classify as locomotion
- **t_thresh = 2 s** - Minimum continuous locomotion duration to filter
- **extra = 30 s** - Additional time after locomotion to exclude (pupil recovery)

### Wheel Parameters

- **Wheel diameter = 15 cm**
- **Encoder resolution = 4028 pulses/rotation**

### Mouse Optical Parameters (Modeling)

- **Focal length = 2.347 mm** (from Remtulla & Hallett)
- **Pupil radius range = 0.3-0.9 mm**
- **Dioptric power = 1/focal_length × n_vitreous**
- **Ocular media transmittance = 0.7**
- **SCE sigma** - Calculated from 12° FWHM (estimated from ground squirrel data)

## Troubleshooting

### Common Issues

**1. "Undefined function or variable 'mouseCore'"**
- This function is not included in the repository
- Contact the authors for access
- Only affects modeling scripts; analysis scripts work without it

**2. "Cannot find isetbio functions"**
- Install isetbio from https://github.com/isetbio/isetbio
- Update the `addpath` line at the top of modeling scripts
- Only affects modeling scripts

**3. "File selection dialog appears empty"**
- Ensure your data files are in MATLAB `.mat` format
- Check that you're in the correct directory
- Verify file permissions

**4. "Array dimensions do not match"**
- Verify your data files match expected format
- Check that frame rates match your experimental setup
- Adjust `stim_fps` and `movie_fps` parameters if needed

### Getting Help

For issues specific to the code:
- Open an issue on GitHub: https://github.com/kerschensteinerd/Fitzpatrick2024_pupillary-contrast-response/issues

For scientific questions:
- Refer to the published manuscript
- Contact the corresponding author

## Citation

If you use this code in your research, please cite:

```bibtex
@article{fitzpatrick2024pupillary,
  title={A pupillary contrast response in mice and humans - Neural mechanisms and visual functions},
  author={Fitzpatrick, et al.},
  journal={[Journal Name]},
  year={2024},
  doi={[DOI if available]}
}
```

Also cite the code repository:

```bibtex
@software{fitzpatrick2024code,
  title={Pupillary Contrast Response Analysis Code},
  author={Fitzpatrick, et al.},
  year={2024},
  url={https://github.com/kerschensteinerd/Fitzpatrick2024_pupillary-contrast-response}
}
```

## Contributing

This is research code accompanying a published manuscript. For bug reports or suggestions, please open an issue on GitHub.

## License

[To be determined - please specify license]

## Acknowledgments

- isetbio toolbox: https://github.com/isetbio/isetbio
- Natural scene data and optical modeling based on established vision science methods

## Contact

For questions about this code or the associated research:
- Open an issue on GitHub
- Contact: [Corresponding author contact information]

---

*Last updated: February 2024*
