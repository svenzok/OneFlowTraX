# OneFlowTraX
OneFlowTraX is a MATLAB app to analyze single-molecule data (e.g. from sptPALM), including localization, tracking, mobility analysis, and cluster analysis. The implementation into one workflow with the visualization of intermediate results facilitates the selection of suitable parameter settings throughout the app. Raw data files can also be batch processed, collecting all parameter settings and results in a comprehensive results file to simplify data management according to FAIR (Findable, Accessible, Interoperable, Reusable) principles. 

# Requirements

For the stand-alone version (v1.3)
- Microsoft Windows 10 or newer, 64-bit
- has to use (older) MATLAB Runtime Version 2025b (R2026a seems to have some .dll dependency issues)
- support for newer graphics cards (e.g. NVIDIA GeForce RTX 5060) only available in R2026a

For execution in MATLAB (v1.4)
- MATLAB R2026a or newer
  - Signal Processing Toolbox 26.1
  - Image Processing Toolbox 26.1
  - Statistics and Machine Learning Toolbox 26.1
  - Parallel Computing Toolbox 26.1 (GPU fitting)

For GPU fitting (Localization step)
  - CUDA capable graphics card with minimum compute capability of 5.2

# How to run
OneFlowTraX v1.3_installer.exe installs the stand-alone version.

OneFlowTraX v1.4.mltbx installs the software as an app in MATLAB (Apps ribbon > MY APPS).
It will actually appear as an add-on. Open the Add-On Explorer, it should appear also here under "Installed" (and can also only be uninstalled from here). Navigate to the OneFlowTraX entry, open the options (three vertical dots), and click "Open Folder". If there are grayed out folders, right-click “Add to Path” > “Selected Folders and Subfolders”.

When copying the entire folder structure and navigating to the main folder (right-click “Add to Path” > “Selected Folders and Subfolders”) the app can be started with OneFlowTraX.mlapp (version 1.4) inside MATLAB's App Designer.

Detailed instructions for installation and the usage of OneFlowTraX can be found in OneFlowTraX v1.0 User Guide.pdf.
You can download sample files and a corresponding noise file from:
https://drive.google.com/drive/folders/12FSUjpxNaUjJI9i2XCgpqXz_B-BYsNJS?usp=sharing

There you will also find the noise file for processing data recorded with EMCCD cameras.

# Copyright and Software License
Copyright (c) SFB1101, ZMPB and IPTC, University of Tübingen, Tübingen.

The scripts of OneFlowTraX are licensed under the [GNU GPL](https://www.gnu.org/licenses/).

# How to cite OneFlowTraX
If you use any of the scripts in OneFlowTraX to process your data, please, cite our [paper](https://doi.org/10.3389/fpls.2024.1358935):
- Leander Rohr, Alexandra Ehinger, Luiselotte Rausch, Nina Glöckner Burmeister, Alfred J. Meixner, Julien Gronnier, Klaus Harter, Birgit Kemmerling, and Sven zur Oven Krockhaus. OneFlowTraX: a user-friendly software for super-resolution analysis of single-molecule dynamics and nanoscale organization. DOI: [https://doi.org/10.3389/fpls.2024.1358935](https://doi.org/10.3389/fpls.2024.1358935)
