# Introduction

In this example a complicated FoV is analyzed using a 1D, 2D and 3D analysis.

# Data Preparation

An example of how to prepare and generate background is found in [PrepAnalysis.ipynb](PrepAnalysis.ipynb). The files generated can be downloaded from [here](https://drive.proton.me/urls/JSGH2TARWW#sgLZ8E2GvsBZ) using the standard VERITAS password.

# 1D Analysis 

A reflected region analysis is performed for [LS 5093 (point source)](./Spectrum_LS5039.ipynb) and [HESS J1825-137 (extended)](./Spectrum_HESS.ipynb).

# 2D Analysis

A [RBM analysis](./RBM_Analysis.ipynb) is performed on the FoV with LS 5093 and HESS J1825-137. An exclusion region is used to mask out the known source and significance maps are obtained for point-like and extended cuts.

# 3D Aanlysis

A [3D Analysis](./3DAnalysis.ipynb) is performed on the entire FoV including LS 5093, HESS J1825-137 and a constant background component. The significances and spectra of the sources are estimated.
