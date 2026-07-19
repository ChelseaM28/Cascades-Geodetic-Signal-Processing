### covariance_realism.py
# This script will quantify the uncertainty of my velocity estimates using a distribution comparison 
# of Monte-Carlo samples of uncertainty realism metrics against the matching chi-squared distribution. 


### Motivating publications:
# Working Group on Covariance Realism. (n.d.). Covariance and Uncertainty 
#        Realism in Space Surveillance and Tracking (A. B. Poore, J. M. Aristoff, 
#        & J. T. Horwood, Eds.). Air Force Space Command Astrodynamics 
#        Innovation Committee.
# Zaidi, Waqar H. , and Matthew D. Hejduk. Earth Observing System Covariance Realism. American Institute 
#        of Aeronautics and Astronautics, 1 Mar. 2016.
# Bos Tero, Machiel Simon, et al. “Introduction to Geodetic Time Series Analysis.” Geodetic Time Series Analysis in 
#        Earth, Springer Verlag, Aug. 2019.
#
#
# Please NOTE that Poore et al quantifies uncertainty against REAL truth, while my uncertainty can
# only be validated against against *simulated* truth. Claiming true uncertainty realism in the context
# of orbit determination, I'd be making an overstatement. However, in the context of GNSS geodesy, 
# my data is adequately 'realistic' since I'm not assuming white noise. TODO: add source.


# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Step 1: Load libraries and data
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
# Step 1: Construct Covariance Matrix from OLS
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -



