# Neurovascular coupling and bilateral connectivity during NREM and REM sleep

This document outlines the steps necessary to generate Kevin L. Turner's data and code for the manuscript **Neurovascular coupling and bilateral connectivity during NREM and REM sleep** by K.L. Turner, K.W. Gheres, E.A. Proctor and P.J. Drew.

---
## Generating the figures
---
This data and code generates all main and supplemental figures and tables that involved data analysis.

Begin by downloading the entire code repository and the data from the following locations:
* Code repository location: https://github.com/DrewLab/Turner_Gheres_Proctor_Drew_Manuscript2020
* Data repository location: https://psu.box.com/s/1yg60jv8ixlzr0ne6ub7v42jexzxnwsu

The box folder contains a pre-analyzed **AnalysisResults.mat** structure that can be used to immediately generate the figures without re-analyzing any data. Download this file (~1.2 GB) as well as the entire github code repository. Add the **AnalysisResults.mat** file to the MATLAB file path by dragging it into the github folder after unzipping the download. Open the MATLAB function **MainScript_Manuscript2020.m** and click the play button to generate the figures. This will take several minutes.

**Software/System Requirements:** Analysis was written in MATLAB 2019b with all optional packages downloaded. Generating the figures and tables are not computationally intensive and will run on most modern laptops and computers. Running **MainScript_Manuscript2020.m** took < 5 minutes to run on a 2018 Macbook Pro (2.6 Ghz 6-Core Intel i7 with 16 Gb 2400 MHz DDR4 RAM, Radeon Pro 560X GPU).

If you would like to automatically save the MATLAB figures and statistical read-outs, change line 29 of **MainScript_Manuscript2020.m** to *saveFigs = 'y';* This will increase the analysis time to create a new folder */Summary Figures and Structures/MATLAB Analysis Figures/* with all of the saved MATLAB figures, PDFs, and text files with key values and statistics. These PDFs were reshaped and stylized in Adobe Illustrator to professionalize each figure.

If you would like to run the core data analysis from the beginning, you will need to obtain all of the original data (~2.5 TB) and maintain the current folder structure on Box. This can be done by downloading all of the data from Box, however, we do not recommend this as Box has limits on the size of downloadable files. Due to these limits on download size (15 GB), this would take hundreds of separate downloads to obtain all the data. We recommend sharing the data through the mailing of a physical copy via external HDD. Please feel free to contact the corresponding author for shipping information and any other questions regarding data sharing. The easiest way to run the analysis from the beginning is to delete or remove the **AnalysisResults.mat** file from the MATLAB file path. If this file is not present, the analysis will attempt to run from the start to recreate it.

Running the core data analysis from the beginning is **highly** computationally expensive and will not run on most computers. We recommend running it on a workstation or analysis computer with at least 64 GB of RAM. **MainScript_Manuscript2020.m** took ~ 48 Hours to run on a 2018 Dell Precision (5820) Workstation (Intel Xeon W-2145 CPU, 64 GB 2666 MHz DDR4 ECC RAM, Nvidia Quadro P1000 GPU)

---
## Original data and pre-processing
---
The data provided on Box has gone through several pre-processing steps. Original data (.TIFF stacks, analog .txt files, camera .bin files, and LabVIEW TDMS files) exceed 17 TB and are available upon request. The code used to initially process all initial files is provided in the code repository's **Pre-processing-scripts** folder. The analysis follows past techniques from Winder et al, 2017. Paper available at: https://www.nature.com/articles/s41593-017-0007-y and code available at https://github.com/DrewLab/Winder_Echagarruga_Zhang_Drew_2017_Code 

Pre-processing is separated into three separate stages for both IOS and 2PLSM data, as well as a series of functions dedicated to sleep scoring for training the machine learning models.

LabVIEW code used to acquire the data can be found at: https://github.com/DrewLab/LabVIEW-DAQ 

---
## Acknowledgements
---
* multiWaitbar.m Author: Ben Tordoff https://www.mathworks.com/matlabcentral/fileexchange/26589-multiwaitbar-label-varargin
* colors.m Author: John Kitchin http://matlab.cheme.cmu.edu/cmu-matlab-package.html
* Chronux subfunctions http://chronux.org/
* Ternary plot subfunctions: Uli Theune https://www.mathworks.com/matlabcentral/fileexchange/7210-ternary-plots
* Several functions utilize varying bits of code written by Dr. Patrick J. Drew and Dr. Aaron T. Winder https://github.com/DrewLab

#### Feel free to contact Patrick Drew or Kevin Turner (klt8@psu.edu) with any issues running the anaysis. 
