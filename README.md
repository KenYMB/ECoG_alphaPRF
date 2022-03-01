# ECoG_alphaPRF
MATLAB scripts and functions from "Spatial Tuning of Alpha Oscillations in Human Visual Cortex"

# Initial Setup  
Please set up toolbox by following the instructions below.
## Use ToolboxToolbox (Recomended)
0. Install [ToolboxToolbox](https://github.com/ToolboxHub/ToolboxToolbox)  
1. Prepare json configuration files by following either instructions below.
   1. Download [Winawerlab toolbox registry](https://github.com/WinawerLab/ToolboxRegistry) for ToolboxToolbox   
      See also: https://wikis.nyu.edu/display/winawerlab/ToolboxToolbox  
      ***Example code for MacOS terminal***
      ~~~
      rm -rf ~/MATLAB/toolboxes/ToolboxRegistry
      git clone https://github.com/WinawerLab/ToolboxRegistry ~/MATLAB/toolboxes/ToolboxRegistry
      ~~~  
   2. Download json configuration files and copy the json files to ToolboxToolbox configurations directory  
      (We will provide download site to get json files only. Currently please get the entire toolbox.)  
      ***Example code for MacOS terminal***
      ~~~
      git clone https://github.com/KenYMB/ECoG_alphaPRF ~/MATLAB/toolboxes/ECoG_alphaPRF
      cp ECoG_alphaPRF/ToolboxToolbox/*.json ~/MATLAB/toolboxes/ToolboxRegistry/configurations/
      ~~~  
2. Solve dependencies  
   ~~~
   tbUse ECoG_alphaPRF
   ~~~  
## Manual Setup (without ToolboxToolbox)
1. Download following Toolboxes and set path to them  
   - ECoG_alphaPRF (This toolbox: https://github.com/KenYMB/ECoG_alphaPRF)
   - ECoG_utils (https://github.com/WinawerLab/ECoG_utils)
   - analyzePRF (https://github.com/cvnlab/analyzePRF)

# Analysis Setup
1. Set working directory  
   Please update `analysisRootPath.m` to set working directory to save data and figures.  
   You can also temporaly update the working directory with the following code.  
   If you don't set 'analysisRootPath', then the toolbox location is considered as the working directory.  
   ~~~
   analysisRootPath = 'ROOT/DIRECTORY/FOR/ANALYSIS';
   ~~~
2. Download data from fileserver  
   (It does not work on Windows)
   ~~~
   s0_downloadData.m
   ~~~
   
# ECoG Analysis
- Sample script `s1_analyzePRF.m` is available.

# Reproduce Figures
- Sample script `s6_makeFigures.m` is available.
