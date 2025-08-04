# Code to produce results from xFuseMap paper
 
# For details see main publication:
A. Hammer et al., ‘Fusion of automatically learned rhythm and morphology features matches diagnostic criteria and enhances AI explainability’, npj Artificial Intelligence, vol. 1, 2025, doi: doi.org/10.1038/s44387-025-00022-w.
 
# Copyright
Alexander Hammer
Institute of Biomedical Engineering
TU Dresden
01307 Dresden, Germany

# Version 1.0, Dresden 31.07.2025
# CC-BY-4.0

#### STRUCTURE ####
## calcRelevanceForBeatIntervals.mat ##
...main script to plot xFuseMaps and do statistics on differences between the distribution of relevance values across intervals in the ECG on beat level 


## \alg ##
...contains scripts that are necessary to run the main script
# .\calculateBoxMetrics.mat
# .\colorbarPlot.mat
# .\defineParameters.mat
# .\doANOVA.mat
# .\highlightSigNiv2.mat
# .\reshapeMultcompareToMatrix.mat
# .\setRelvToNan.mat
# .\sortVars.mat


## \data ## 
\data\xFuseMap_test_full.mat
... contains ECGs and relevance information of the xFuseMap® test set, containing n=1521 ECGs, classified using xECGArch® and explained using deep Taylor decomposition
# db 		      database information for each recording (1 x n string)
# fp_labels 	fiducial points labels (1 x 9 string)
# fp_locs 	  fiducial points locations for each recording (n x 1 cell)
# fw_locs 	  f wave locations for each recording (n x 1 cell)
# int_labels	interval labels (n x 10 string)
# int_locs 	  interval locations (n x 1 cell)
# labels 		  label annotation [1: AF, 0: n-AF] (1 x n double)
# labels_full label annotatio [1: AF, 2: NSR, 0: O) (1 x n double)
# meta 		    meta data (age, sex, ...) (n x 5 table)
# pred 		    model prediction certainty (4 x n double)
# pred_weight model prediction certainty weighted (4 x n double)
# rec 		    recording ID (1 x n string)
# relv 		    relevance values (5000 x n x 4 double)
# relv_int 	  relvance per interval (1 x 1 struct)
# relv_meth 	relevance extraction method (xAI="dt") (1 x 1 string)
# sigs 		    ECG recordings (preprocessed) (5000 x n double)
# vn_pred 	  variable names: model prediction (4 x 1 string)
# vn_relv  	  variable names: relevance (4 x 1 string)


## \results ##
...folder in which the xFuseMaps and statistical results are saved when running the main script.
