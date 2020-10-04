# Neurovascular coupling and bilateral connectivity during NREM and REM sleep

This document outlines the steps necessary to generate the figures for the manuscript **Neurovascular coupling and bilateral connectivity during NREM and REM sleep** by K.L. Turner, K.W. Gheres, E.A. Proctor and P.J. Drew.

---
## Generating the figures
---
This data and code generates all main and supplemental figures and tables that involved data analysis.

Begin by downloading the entire code repository and the data from the following locations:
* Code repository location: https://github.com/DrewLab/Turner_Gheres_Proctor_Drew_eLife2020
* Data repository location: https://doi.org/10.5061/dryad.6hdr7sqz5     

The Dryad link contains a pre-analyzed **AnalysisResults.mat** structure that can be used to immediately generate the figures without re-analyzing any data. Download this file (~1.2 GB) as well as the entire github code repository. Add the **AnalysisResults.mat** file to the MATLAB file path by dragging it into the folder containing the code. Open the MATLAB function **MainScript_eLife2020.m** and run.

**Software/System Requirements:** Code was tested with MATLAB 2019b. Running **MainScript_eLife2020.m** took < 5 minutes to run on a 2018 Macbook Pro (2.6 Ghz 6-Core Intel i7 with 16 Gb 2400 MHz DDR4 RAM, Radeon Pro 560X GPU).

If you would like to automatically save the MATLAB figures and statistical read-outs, change line 29 of **MainScript_eLife2020.m** to *saveFigs = 'y';* This will increase the analysis time to create a new folder */Summary Figures and Structures/MATLAB Analysis Figures/*.

LabVIEW code used to acquire the data can be found at: https://github.com/DrewLab/LabVIEW-DAQ

---
## Acknowledgements
---
* multiWaitbar.m Author: Ben Tordoff https://www.mathworks.com/matlabcentral/fileexchange/26589-multiwaitbar-label-varargin
* colors.m Author: John Kitchin http://matlab.cheme.cmu.edu/cmu-matlab-package.html
* Chronux subfunctions http://chronux.org/
* Ternary plot subfunctions: Uli Theune https://www.mathworks.com/matlabcentral/fileexchange/7210-ternary-plots
* Several functions utilize code written by Dr. Patrick J. Drew and Dr. Aaron T. Winder https://github.com/DrewLab

#### Contact Patrick Drew (pjd17 psu edu) or Kevin Turner (klt8 psu edu) for further information.
