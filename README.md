# ECoG_alphaPRF

Analysis tools and MATLAB scripts for Yuasa et al., 2023, "Spatial tuning of alpha oscillations in human visual cortex" (eLife reviewed preprint: https://elifesciences.org/reviewed-preprints/90387). This repository contains all necessary code to reproduce the analyses and figures from the paper.

# Initial Setup  

## Option 1: Using ToolboxToolbox (Recommended)

1. Install [ToolboxToolbox](https://github.com/ToolboxHub/ToolboxToolbox)

2. Set up the configuration files using either method:

   ### Method A: Using Winawerlab Registry
   1. Download the [Winawerlab toolbox registry](https://github.com/WinawerLab/ToolboxRegistry):
      ```bash
      # For MacOS
      rm -rf ~/MATLAB/toolboxes/ToolboxRegistry
      git clone https://github.com/WinawerLab/ToolboxRegistry ~/MATLAB/toolboxes/ToolboxRegistry
      ```
      See also: https://wikis.nyu.edu/display/winawerlab/ToolboxToolbox

   ### Method B: Direct Configuration
   1. Download and copy configuration files:
      ```bash
      # For MacOS
      git clone https://github.com/KenYMB/ECoG_alphaPRF ~/MATLAB/toolboxes/ECoG_alphaPRF
      cp ECoG_alphaPRF/ToolboxToolbox/*.json ~/MATLAB/toolboxes/ToolboxRegistry/configurations/
      ```

3. Install dependencies:
   ```matlab
   tbUse ECoG_alphaPRF
   ```

## Option 2: Manual Setup

1. Download and set path to the following toolboxes:
   - ECoG_alphaPRF (This toolbox: https://github.com/KenYMB/ECoG_alphaPRF)
   - ECoG_utils (https://github.com/WinawerLab/ECoG_utils)
   - analyzePRF (https://github.com/cvnlab/analyzePRF)
   - Fieldtrip (https://github.com/fieldtrip/fieldtrip)

# Analysis Setup

## 1. Set Working Directory

You have three options to set the working directory for saving data and figures:

A. Edit the file `analysisRootPath.m` in the 'Environment' folder:
   ```
   rootPath = 'YOUR/ANALYSIS/DIRECTORY';  % Edit this line in the file
   ```

B. Set the analysis directory by running this command in MATLAB:
   ```matlab
   analysisRootPath = 'YOUR/ANALYSIS/DIRECTORY';  % The toolbox will use this directory for analysis
   ```
   Note: This sets a MATLAB variable that the toolbox uses to determine where to save data and figures.

C. Copy `analysisRootPath.m` to your current working directory (the toolbox will use this location)

Note: If no working directory is specified through any of these methods, the toolbox defaults to its installation location.

## 2. Download Data

The BIDS-formatted dataset is available on OpenNeuro: [Visual ECoG dataset](https://openneuro.org/datasets/ds004194)

To automatically download the data (not available on Windows):
```matlab
s0_downloadData
```

# Getting Started

- Example scripts are available in the 'SampleScripts' directory
- To reproduce the paper's main figures, run `s7_makeFigures`
- For supplementary figures, run `s8_makeSupplementaryFigures`
