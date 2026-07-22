### covariance_realism.py

# !! ----------- WORK IN PROGRESS ----------- !!


# This script will quantify the uncertainty of my velocity estimates using a distribution comparison 
# of Monte-Carlo samples of uncertainty realism metrics (at n = 1) against the matching chi-squared distribution. 
# Jul 18, 2026
# Version 1
# Chelsea Momoh

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
import math
from scipy.signal import periodogram

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

station_lengths = {
    "p349": p349.nrows(),
    "p380": p380.nrows(),
    "p434": p434.nrows(), 
    "p441": p441.nrows()
}

#TODO: All k values are to be replaced with values from "alphas" list.
station_slopes= {#TODO: All k values are to be replaced with values from "alphas" list.
    "p349_east": k, "p349_north": k, "p349_vert": k,
    "p380_north": k, "p380_east": k, "p380_vert": k,
    "p434_north": k, "p434_east": k, "p434_vert": k,
    "p441_north": k, "p441_east": k, "p441_vert": k,
}
print("Finished loading data")



# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Mathematical Reasoning for Uncertainty Realism Script
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -

#Among everything, be sure to explain why z² is the Mahalanobis Distance Metric from Poore et al. at n = 1

# STEP 1: PARAMETRIC ESTIMATION OF C for GLS


# STEP 3: CREATING COLORED MONTE CARLO NOISE
#To find the σ^2_colored and σ^2_white, I look at the power of each PSD graph at the 
# intercept (where sampling frequency = frequency) and at the flattening point (line of demarcation which
# I already set to ~5 cycles/yr)
#looking back at my PSD graphs, var_white will eb the power located at the frewuency of 10^(0.7)
#var_colored will be the power at the frequency of 365.25 #TODO: MIGHT NOT WORK - aliasing means i dont have this.
#However, this will need to be calculated for each direction for each station.



# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Step 1: Parametric Estimation of C for GLS
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -

#Brief: This section will use define a noise model using discoveries from OLS
def create_general_power_law_cov_matrix(station, k):
    h = [1] 
    n = station_lengths[station]
    for i in range(n): #Make sure there is no mismatch between list beginnign with 0/1
        h.append((i - (k/2) - 1)*(h[i-1]/i)) #This is how h is defined in Tero et al
    H = np.zeros((n,n)) #I am going to convert h into matrix form, H
    for i in range(n): #I need to fix this. it's not updating any values.
        H[k, :k+1] = h[k::-1] #This is creating a lower triangular matrix
    print(f"Completed general power-law covariance matrix for {station}.")
    J = H@np.transpose(h) 
    return H, h, J
#J is a covariance matrix created from colored noise, dependent upon k.
#       when multiplied by var_colored it is the colored component of covariance. 
#       J is not needed for monte carlo.
        


def find_characterized_var(station):
    white_freq = 10**(0.7)
    colored_freq = 365.25 #TODO: MIGHT NOT WORK - aliasing means i dont have this.
    var_white = #i need a function that takes a freq and outputs its power form a FITTED 
    #line of the PSD graph. 

def create_parametric_C(station):
    k = station_slopes[station]
    n = station_lengths[station]
    H, h, J = create_general_power_law_cov_matrix(station, k)
    var_white, var_colored = find_characterized_var(station)
    I = np.eye(n)
    C = var_colored@J + var_white@I
    return C #This C will be used to contruct the GLS equation.


# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Step 2: Describe OLS and GLS equations.
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -

# position(t) = a + b·t + c·sin(2πt) + d·cos(2πt) + e·sin(4πt) + f·cos(4πt) + residual

#@Brief: Plain OLS and GLS Equations (not specific, only general representatives):
#OLS: 
#GLS: 
#The difference between the two is the weighted component of GLS and its assumptions, etc etc...

#@Brief: This section will verify all assumptions for my methodology are met before running Monte Carlo.


# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Step 3: Use Monte Carlo to Generate N Synthetic Ground Station Series
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -

#STEP by STEP: 
#       I already have C, which is set permanently. 
#       Hv = w #This is the filter used to create a vector of colored noise for a simulation.
#.      w is the vector of noise we need for each simulation. 
#       H = summation(h_i) written in matrix format.
#       h_i = binom(n,k) function and it is a set value for a simulation.
#       v is a vector with independent and identically distributed (IID) Gaussian noise
#       
#       

#The goal of this section is to generate N Synthetic ERROR series
def generate_monte_carlo_series(H, station):
    directions = ["north", "east", "vert"]
    noise_models = []
    all_synthetic_series = []

    for direction in directions:
        station_direction = str(station) + "_" + str(direction)
        #This sets the TRUE velocity that OLS and GLS will conform to
        VELOCITY = betas[station_direction][1].round() #I need to think about whether this is deterministic
        a, c, d, e, f = betas[station_direction][[0, 2, 3, 4, 5]]
        N = 500 #subject to change.
        size = station_lengths[station]
        
        #Brief: This section generate N synthetic ground station motion time series.
        for simulation in range(N):
            v= np.random.normal(loc=0.0, scale=1.0, size=size)
            w = H@v #summation(h)@v
            noise_models.append(w) #I will only begin saving persistently once I confirm the code works!
        #synthetic_series = [signal] + [noise]
        #synthetic_series = [X][B_true] + [noise]
            synthetic_series = X_matrices@[a, VELOCITY, c, d, e, f] + w #TODO: ALL of the X matrices? No. Not all of them. Fix this.
            all_synthetic_series.append(synthetic_series)

    return all_synthetic_series

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


