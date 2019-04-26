# riksneurotools
Some Matlab functions from Rik Henson. Use at your own risk!

GLM contains functions related to the General Linear Model (particularly for fMRI)

SPM contains functions specific to SPM (eg bypassing 2nd-level batch interface)

Conn contains functions for functional connectivity (particularly for fMRI)

Util contains utility functions

MEEG contains functions for MEEG analysis (in SPM)

GLM
===
GLM directory contains some functions for General Linear Model on single vector of data (glm.m), a repeated-measures ANOVA (repanova.m), independent and repeated measures T-tests (t_matrix), error nonsphericity in mass univariate ANOVAs (check_pooled_error.m), efficiency of fMRI GLM (fMRI_GLM_efficiency.m), GLM for resting-state fMRI functional connectivity (rsfMRI_GLM.m) and GLMs for estimating single fMRI trials using LSA/LSS (fMRI_multitrial_GLMs.m).

SPM
===
SPM directory contains functions for SPM (only one at moment for generic ANOVAs, bypassing GUI)

Conn
====
rsfMRI_GLM.m estimates connectivity for fMRI, including covariates, filtering and prewhitening
connectivity_stats.m just tests connectivity matrices between groups, with various thresholding schemes
dcor_dc.m and dcor_uc.m estimate distance correlation between to two multivariate measures (eg two ROIs)

Util
====
roi_extract.m extracts mean (or first singular vector) across voxels within ROIs defined by mask image;

MEEG
====
MEEG directory contains functions for M/EEG analysis in SPM (only one at moment for ICA-cleaning)

