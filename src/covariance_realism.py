### covariance_realism.py
# This script will quantify the uncertainty of my velocity estimates using a distribution comparison 
# of Monte-Carlo samples of uncertainty realism metrics (at n = 1) against the matching chi-squared distribution. 


### Motivating publications:
# Bos Tero, Machiel Simon, et al. “Introduction to Geodetic Time Series Analysis.” Geodetic Time Series Analysis in 
#        Earth, Springer Verlag, Aug. 2019.
# Working Group on Covariance Realism. (n.d.). Covariance and Uncertainty 
#        Realism in Space Surveillance and Tracking (A. B. Poore, J. M. Aristoff, 
#        & J. T. Horwood, Eds.). Air Force Space Command Astrodynamics 
#        Innovation Committee.
# Zaidi, Waqar H. , and Matthew D. Hejduk. Earth Observing System Covariance Realism. American Institute 
#        of Aeronautics and Astronautics, 1 Mar. 2016.
#
#
# Please NOTE that Poore et al quantifies uncertainty against REAL truth, while my uncertainty can
# only be validated against against *simulated* truth. Claiming true uncertainty realism in the context
# of orbit determination, I'd be making an overstatement. However, in the context of GNSS geodesy, 
# my data is adequately 'realistic' since I'm not assuming white noise. TODO: add source.


# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Step 0: Load libraries and data
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
import os
import time
os.chdir("/workspaces/GNSS/data")
import pandas as pd
import numpy as np  
import json
import matplotlib.pyplot as plt
#For PELT change point detection
import ruptures as rpt

with open("alphas.json", "r") as f:
    alphas = json.load(f)
with open("betas.json", "r") as f:
    betas = json.load(f)
betas = {key: np.array(value) for key, value in betas.items()}
with open("residuals.json", "r") as f:
    residuals = json.load(f)
residuals = {key: np.array(value) for key, value in residuals.items()}
with open("X_matrices.json", "r") as f:
    X_matrices = json.load(f)
X_matrices = {key: np.array(value) for key, value in X_matrices.items()}
with open("metadata.json", "r") as f:
    metadata = json.load(f)
with open("changepoints.json", "r") as f:
    changepoints = json.load(f)
changepoints = {key: np.array(value) for key, value in changepoints.items()}


p349 = pd.read_json("p349.json", orient="records")
p349['Date'] = pd.to_datetime(p349['Date'])
p380 = pd.read_json("p380.json", orient="records")
p380['Date'] = pd.to_datetime(p380['Date'])
p434 = pd.read_json("p434.json", orient="records")
p434['Date'] = pd.to_datetime(p434['Date'])
p441 = pd.read_json("p441.json", orient="records")
p441['Date'] = pd.to_datetime(p441['Date'])

stations = {'p349': p349, 'p380': p380, 'p434': p434, 'p441': p441}
station_dates = {
    "p349": p349['Date'], "p380": p380['Date'], "p434": p434['Date'], "p441": p441['Date'],
}
print("Finished loading data")



# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Mathematical Reasoning for Uncertainty Realism Script
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -

#Among everything, be sure to explain why z² is the Mahalanobis Distance Metric from Poore et al. at n = 1

# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Step 1: Parametric Estimation of C for GLS
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -

#Brief: This section will use define a noise model using discoveries from OLS


# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Step 2: Define OLS and GLS equations.
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -

#@Brief: Plain OLS and GLS Equations (not specific, only general representatives):
#OLS: 
#GLS: 
#The difference between the two is the weighted component of GLS and its assumptions, etc etc...

#@Brief: This section will verify all assumptions for my methodology are met before running Monte Carlo.


# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Step 3: Use Monte Carlo to Generate N Synthetic Ground Station Series
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -

#Brief: This section will set the TRUE velocity that OLS and GLS will conform to
VELOCITY = _
N = _

#Brief: This section generate N synthetic ground station motion time series.

# * - * - * - * - * - * - * - * - * - * - * - *
# Step 4: Refit OLS/GLS on each Series   
# * - * - * - * - * - * - * - * - * - * - * _ *

#@Brief: Data manipulation to format for matrix calculations

#@Brief: This section will recover velocity estimates for OLS and GLS from each simulation

#Brief: This section will recover σ²_OLS and σ²_GLS from each simulation


# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Step 5: Create Normalized σ² Vectors (Uncertainty Realism Metrics) for Comparison 
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# It is helpful to recall that z² is the Mahalanobis Distance Metric from Poore et al. at n = 1


# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Step 6: Plot Velocity Estimates on Histogram (Not necessary, purely as visual additions)
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -


# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Step 7: Distribution Comparison (FOCAL POINT OF PROJECT)
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -

#@Brief: This section will plot OLS/GLS z² values on a histogram


#@Brief: This section will perform a goodness of fit test to confirm visual results


# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Analysis of the significance of the OLS/GLS σ² histograms
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -




# What a doozy!

# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Future Modifications For this project
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -

# 1. Array Procesing Extension across all stations
# 2. Utilize MLE rather than periodogram slope to define alpha 
#       (also k in some literature) to build parametric C. See Tero et al.
# 3. Pull USGS earthquake catalogue to confirm changepoints
# 4. Segment residual regimes between changepoints, then flag outliers (IQR/ Z-score)
# 5. Model diagnostics (Model-order, multicollinearity, etc)
# 6. Wrap project as real-time dynamic pipeline


