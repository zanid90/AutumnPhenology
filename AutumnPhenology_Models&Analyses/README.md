# This folder contains all scripts of correlation analyses, past and future autumn phenology predictions, model performance analyses, and related output

### The two subfolders contain the following documents and scripts:
- Data:
	- ModelAnalysis_1_Predicted_DoYoff.csv: predictions of leaf senescence dates for 10 process-based autumn phenology models.
	- ModelAnalysis_2_OptimalParameters.csv: optimal parameters at the timeseries-level for the 10 models.
	- ModelAnalysis_3_Performance_R2_RMSE_slope.csv: performance statistics of observed vs. predicted values of autumn anomaly for the 10 models.
	- ModelAnalysis_4_CrossValidation_RMSE.csv: validation statistics of observed vs. predicted values of autumn anomaly for the 10 models.
	- ModelAnalysis_5_FutureAutumnProjections.csv: future predictions of leaf senescence dates for 6 process-based autumn phenology models.
- Codes:
	- ModelAnalyses_1_Correlation&PathAnalyses_Figures1.r: analysing correlations between autumn phenology drivers and leaf senescence dates, and generating related figures (outputs: Fig. 1ABC, Fig. S1, Fig. S2).
	- ModelAnalyses_2_AutumnPhenology_Optimization_Predictions.r: calibrating parameters for the 10 models at the timeseries-level and estimating leaf senescence dates for 1948-2015 (outputs: ModelAnalysis_1_Predicted_DoYoff.csv and ModelAnalysis_2_OptimalParameters.csv).
	- ModelAnalyses_3_Models_Performance_Figures2.r: evaluating the 10 process-based autumn phenology models, validating the photosynthesis model, and generating related figures (outputs: ModelAnalysis_3_Performance_R2_RMSE_slope.csv, ModelAnalysis_4_CrossValidation_RMSE.csv, Fig. 3ABCD, Fig. S8).
	- ModelAnalyses_4_AutumnPhenology_FutureProjections_Figures3.r: predicting leaf senescence dates for 2016-2100, and generating related figures (outputs: ModelAnalysis_5_FutureAutumnProjections.csv, Fig. 4ABCD, Fig. S6, Fig. S7).
