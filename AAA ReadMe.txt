This is the repository needed to reproduce the results provided in the manuscript
"Evaluating quantile treatment effects with machine learning: An application to the informal sector wage gap".
More specifically, the file "DataPrep.R" allows the reader to open and recode the variables of interest.
Once the data are prepared, the proposed classification models and the quantile regression model can be implemented 
by running the code reported in the file "RunModels.R" while the approximation of the standard errors
of the quantile regression model can be computed by running the code provided within the file "SeQuantile.R". Additionally, the result presented in Figure 1 of the paper can be reproduced by runnig the script contained in the file PlotRes.R
Finally, the simulation study can be reproduced by running the code provided in the file "SimGam.R"
