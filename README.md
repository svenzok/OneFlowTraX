# OneFlowTraX
OneFlowTraX is a MATLAB app to analyze single-molecule data (e.g. from sptPALM), including localization, tracking, mobility analysis, and cluster analysis. The implementation into one workflow with the visualization of intermediate results facilitates the selection of suitable parameter settings throughout the app. Raw data files can also be batch processed, collecting all parameter settings and results in a comprehensive results file to simplify data management according to FAIR (Findable, Accessible, Interoperable, Reusable) principles. 

# Requirements

For the stand-alone version
- Microsoft Windows 10 or newer, 64-bit

For execution in MATLAB
- MATLAB R2023a or newer
  - Signal Processing Toolbox 9.2
  - Image Processing Toolbox 11.7
  - Statistics and Machine Learning Toolbox 12.5
  - MATLAB Parallel Computing Toolbox 7.8 (GPU fitting)

For GPU fitting (Localization step)
  - CUDA capable graphics card with minimum compute capability of 5.2

# How to run
OneFlowTraX_Installer.exe installs the stand-alone version.

OneFlowTraX.mlappinstall installs the software as an app in MATLAB (Apps ribbon > MY APPS).

When copying the entire folder structure and navigating to the main folder (right-click “Add to Path” > “Selected Folders and Subfolders”) the app can be started with OneFlowTraX.mlapp.

Detailed instructions for installation and the usage of OneFlowTraX can be found in OneFlowTraX v1.0 User Guide.pdf.
You can download sample files and a corresponding noise file from:
https://drive.google.com/drive/folders/12FSUjpxNaUjJI9i2XCgpqXz_B-BYsNJS?usp=sharing

# Copyright and Software License
Copyright (c) SFB1101, ZMPB and IPTC, University of Tübingen, Tübingen.

The scripts of OneFlowTraX are licensed under the [GNU GPL](https://www.gnu.org/licenses/).

# How to cite OneFlowTraX
If you use any of the scripts in OneFlowTraX to process your data, please, cite our [paper](https://www.biorxiv.org/content/10.1101/2023.08.10.552827v1):
- Leander Rohr, Alexandra Ehinger, Luiselotte Rausch, Nina Glöckner Burmeister, Alfred J. Meixner, Julien Gronnier, Klaus Harter, Birgit Kemmerling, and Sven zur Oven Krockhaus. OneFlowTraX: A User-Friendly Software for Super-Resolution Analysis of Single Molecule Dynamics and Nanoscale Organization. bioRxiv preprint DOI: [https://doi.org/10.1101/2023.08.10.552827](https://doi.org/10.1101/2023.08.10.552827)
