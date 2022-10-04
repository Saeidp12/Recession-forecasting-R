# Recession Forecasting Using Bayesian Classification
# Paper implementation in R

(This is a repository for one of my freelance projects)

The paper and its supplemantary MATLAB code were downloaded from the official website for Federal Reserve Bank of Kansas City.

In this repository, the matlab code files from a research paper are converted to R. The paper is about a time-series analysis.

Step 1: data is transformed based on the transform code. 

Step 2: Values far above or far below median are set to NA (Not Available: missing).

Step 3: Data is lagged for 10 months. 

Step 4: It doesn’t use the last 18 months in the forecast (It is blind to them: forecast_blind).

Step 5: It selects dates that are common between macro data and recession data. 

Step 6: The starting date for forecast is 1973/01/01. 

Step 7: The ending date for forecast is the last common date between macro and recession which is 2017-01-01.

Step 8: For every date between 1973/01/01 and 2017/01/01 this process happens for forecasting:

    For example, if the date is 2010/01/01, all the data from the beginning of the common data frame (1960/01/01 up to (2010/01/01 – 1m – 18m – forecast_horizon) is used as train data. If the horizon window is 3 months (22 months total) it should be up to 2008/01/03. Thus, for this train data which starts from 1960/01/01 to 2008/01/03, the forecast should be for the date 2010/01/01. 

In the forecast table, the element in row with the date 2010/01/01 and column H3 belongs to this forecast. 
