spatialNetSmooth_scripts
=========
This repository contains the scripts for reproducing the plots from my thesis.
# Datasets
The Visium breast cancer sample can be found on the 10xGenomics web-
site: https://www.10xgenomics.com/datasets/. The sample used was Human Breast Cancer (Block A Section 1).
The file "tissue_lowres_image_annotated.png" should be downloaded from this repo.
The second dataset belongs to a paper by Anderson et al. (2021) and
also derives from breast cancer samples. The processed count matrices, HE-
images and metadata with spot-annotation can be found here: https://zenodo.
org/records/4751624. The samples used in this thesis were A1, B1, C1, D1,
E1 and H1.


# How to use:
### Reproduce Spatial plots for Visium dataset:
Run "run_tests.R" and specify input and output directory. After that run "calc_R2_spatialplots_visium.R" and look inside the "_R2.csv"-files to find the best smoothing parameters for each method. Then, the column from the score dataframe can be selected in the script and used for plotting.

### Reproduce Spatial plots for ST dataset: 
Run "run_tests_ST.R" for each sample and change "letter" to the current sample. Specify input and output directory. In the output directory, there should be folders named like the samples (e.g. "A1").
After that run "calc_R2_spatialplots_ST.R" (also for every sample) and look inside the "_R2.csv"-files to find the best smoothing parameters for each method. Then, the column from the score dataframe can be selected in the script and used for plotting.

### Reproduce plot of F1-scores for Visium data:
Run "run_tests.R" and specify input and output directory. Then look into the F1-plots for every method and choose the best smoothing paramters. Specify the column numbers at the bottom of the script for plotting "F1_best".

### Reproduce table comparing accuracy of thresholding methods for Visium dataset:
Run "run_test.R" and then "script_accuracy_calc_Visium.R"

### Reproduce R²-compare table/plot:
Run "calc_R2_spatialplots_visium.R" and "calc_R2_spatialplots_ST.R" (for each sample) and look inside the "_R2.csv"-files to find the best R²-score for each method per sample. Then run "R2_plot.R".
