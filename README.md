# riksneurotools
Some Matlab functions from Rik Henson.

roi_extract.m extracts mean (or first singular vector) across voxels within ROIs defined by mask image

GLM directory contains some functions for General Linear Model on single vector of data (glm.m), a repeated-measures ANOVA (repanova.m), independent and repeated measures T-tests (t_matrix), error nonsphericity in mass univariate ANOVAs (check_pooled_error.m), efficiency of fMRI GLM (fMRI_GLM_efficiency.m), GLM for resting-state fMRI functional connectivity (rsfMRI_GLM.m).

SPM directory contains functions for SPM (only one at moment for generic ANOVAs, bypassing GUI)

MEEG directory contains functions for M/EEG analysis in SPM (only one at moment for ICA-cleaning)

