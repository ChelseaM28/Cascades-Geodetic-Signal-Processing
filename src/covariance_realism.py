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
#
# Please NOTE I must undergo a test to confirm the true GNSS stations' residuals are gaussian to align with the
# assumptions maded in Poore et al. Otherwise I would need to generalize from a covariance realism metric (which I 
# initially planned to do) to an uncertainty realism metric, which relaxes the gaussian assumption. Though my Monte
# carlo simulation will remain internally consistent whichever route I choose, recall, as stated in the note above,
# simulated truth is not necessarily real truth. In other words, for my results to be relevant, the assumption must 
# match reality (which is kinda the core theme of this whole project!!)


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
    "p349": len(p349),
    "p380": len(p380),
    "p434": len(p434), 
    "p441": len(p441)
}

#NOTE THESE ARE NOT ACCURATE YET. They must be calculated using the fitted line i will use for the find_characterized_var function
# I should also probably segment according to the changepoints I calculated earlier. Means the alphas list I generated earlier will likely be scratched.
station_slopes= {
    "p349_north": alphas["p349_north"], "p349_east": alphas["p349_east"], "p349_vert": alphas["p349_vert"],
    "p380_north": alphas["p380_north"], "p380_east": alphas["p380_east"], "p380_vert": alphas["p380_vert"],
    "p434_north": alphas["p434_north"], "p434_east": alphas["p434_east"], "p434_vert": alphas["p434_vert"],
    "p441_north": alphas["p441_north"], "p441_east": alphas["p441_east"], "p441_vert": alphas["p441_vert"],
}
print("Finished loading data")



# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Mathematical Reasoning for Uncertainty Realism Script
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -

# MAIN PREMISE: We are determining whether OLS or GLS has more accurate formal uncertainty (σ^2). In order to do this, we need a
# GROUND TRUTH to verify against (in my case, ground truth is simulated). To determine whether σ^2_OLS or σ^2_GLS is better, 
# use Monte Carlo to generate N synthetic series with a predetermined VELOCITY, and some synthetic noise (all of which 
# exhibiting the same noise model but with different INDIVIDUAL noise values).
# Don't be confused: We refit OLS and GLS on each Monte Carlo series. When we recover the model's estimated B values from the simulations,
# and plot them on ahistogram, both models will center on VELOCITY_true, but the spreads will be different. That's not the end though.

# Each model should be aware of its own spread (spread ~ uncertainty). So we calculate σ^2 for each simulation, normalize them to a z distribution 
# to make them comparable, and plot them on a histogram. 
# (THESE normalized values, z, are my 'uncertainty realism' metrics). If the σ^2 values are accurate,
# the z distribution should follow a chi-sqrd distribution with df = 1.
#     WHY? - In order to normalize the errors, we are doing z = (B_estimated - B_true)/σ_formal so the Var(z) = σ^2_true/σ^2_formal
#            In other words, we are dividing actual uncertainty by what the model thinks its uncertainty is. Z^2 is actually the Mahalanobis distance at n = 1.
#            This normalization will result in a chi-sqrd distribution. (Mahalanobis distance, btw, is essentially the distance between a point and a distribution. 
#            I wont get in the weeds though.)TODO: STILL NEEDS EXTRA CLARIFICATION on WHY CHI SQRD?
#TODO Among everything, be sure to explain why z² is the Mahalanobis Distance Metric from Poore et al. at n = 1
# The issue we will encounter is OLS is BLIND to colored noise, 
#     WHAT? - In other words, rather than considering correlated errors (σ^2_OLS@C), OLS assume the error is independent, and thus operates 
#             under the premise that σ^2_OLS@C = σ^2_OLS@I. <-- identity matrix 
# so its z distribution will be shifted right/heavier on the tail. This means its own idea of how uncertain it is is wrong. 
#
# As far as I'm aware, theoretically, aside from model uncertainty, with large enough N, the methodologies were built s.t. the standardized errors 
# would always follow the appropriate distribution. But we will see that OLS breaks and any inferences we make with using its reported uncertainty
# (any hypothesis tests or confidence intervals) will NOT be credible.


# SCRIPT STRUCTURE: I divide this script into distinct steps with subsections. A general overview of each step is found below.

# STEP 1: PARAMETRIC ESTIMATION OF C for GLS
# When fitting a model (estimating the parameters, B) generalized least squares utilizes
# a covariance matrix to correct against heteroscedasticity and correlated error.
# My monte carlo simulation needs to generate N synthetic series with identical covariance 
# but different values in each epoch. In this way OLS and GLS can be refit under the same error
# conditions. So once the covariance of the errors (not parameters TODO:CHECK) is determined, 
# it will be reused without any changes. To estimate the covariance accurately, I use Tero et al's
# parametric method of calculating C.


# The parametric method of estimating C requires the colored portion of the variance σ^2_colored be
# multiplied by the appropriate covariance model/matrix, J, to for the colored component of C. e.g. σ^2_colored@J
# Next, the white noise variance, σ^2_white, will be multiplied by its covariance model (which will be an identity matrix).
# e.g. σ^2_white@I
# So C = σ^2_colored@J + σ^2_white@I
# To find the σ^2_colored and σ^2_white, I look at the power of each PSD graph at the 
# intercept (where sampling frequency = frequency) and at the flattening point (line of demarcation which
# I already set to ~5 cycles/yr)
#looking back at my PSD graphs, var_white will eb the power located at the frewuency of 10^(0.7)
#var_colored will be the (extrapolated) power at the frequency of 365.25 
#However, this will need to be calculated for each direction for each station.

# STEP 3/4: CREATING COLORED MONTE CARLO SERIES + REFITTING
# I create N error series (w) 
# and I add them to N deterministic signals (a + t + sin(2πt) + cos(2πt) + sin(4πt) + cos(4πt)))
# which creates the N series OLS and GLS will be refitted to. 

# STEPS 5-6: These will complete the distribution comparison and interpret the results as 
# described above.


# END MATHEMATICAL REASONING



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
        H[i, :i+1] = h[i::-1] #This is creating a lower triangular matrix
    print(f"Completed general power-law covariance matrix for {station}.")
    J = H@np.transpose(np.array(h[:n]))
    return H, h, J
#J is a covariance matrix created from colored noise, dependent upon k. In my text, k seems 
#       synonymous with k. J, when multiplied by var_colored (σ^2_colored), is the colored 
#       component of covariance. Once used to calculate C, J is not needed again for monte carlo.
        

def fit_psd_line(residual, fs=365.25, freq_cutoff=5):
    #o compute alphas.json, on the same restricted range
    #(aka masked frequency range). 
    freqs, power = periodogram(residual, fs=fs)
    bin_centers, bin_means = bin_psd(freqs, power)
    mask = (bin_centers < freq_cutoff) & (~np.isnan(bin_means))
    log_f = np.log10(bin_centers[mask])
    log_p = np.log10(bin_means[mask])
    slope, intercept = np.polyfit(log_f, log_p, 1)
    return slope, intercept


def power_at_freq(freq, slope, intercept):
    #Reads power off the fitted log-log line at a given frequency.
    log_power = slope * np.log10(freq) + intercept
    return 10 ** log_power


def find_characterized_var(station, direction):
    white_freq = 10**(0.7)
    colored_freq = 365.25 #This is the value at which f/f_s = 1, so when plugged 
    # into the equation for power-law noise (figure 22 in Tero et al.), the resulting 
    # value is P, a constant, the power, the variance. Yes, this frequency, f, is beyond 
    # the nyquist frequency. However, if my understanding of the literature is correct, 
    # f always is. However, I needed the aforementioned ratio to hold true. 
    # It's a real idiosyncracy I'm wrestling with.  

    #To fing var_white and var_colored, i needed a function 
    # that takes a freq and outputs its power form a FITTED line of the PSD graph. 
    key = f"{station}_{direction}"
    slope, intercept = fit_psd_line(residuals[key])
    var_white = power_at_freq(white_freq, slope, intercept)
    var_colored = power_at_freq(colored_freq, slope, intercept)
    return var_white, var_colored

def create_parametric_C(station, direction):
    key = f"{station}_{direction}"
    k = station_slopes[key]
    n = station_lengths[station]
    H, h, J = create_general_power_law_cov_matrix(station, k)
    var_white, var_colored = find_characterized_var(station, direction)
    I = np.eye(n)
    C = var_colored*J + var_white*I
    return C #This C will be used to contruct the GLS equation.


# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Step 2: Describe OLS and GLS equations.
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -

# The data is modeled using the following predictors:
# position(t) = a + b·t + c·sin(2πt) + d·cos(2πt) + e·sin(4πt) + f·cos(4πt) + residual

#@Brief: Plain OLS and GLS Equations (not specific, only general representatives):
#OLS: y = XB + e  <-- Very familiar!
#     B = (X'X)^(-1)  <-- Algebraic manipulation to isolate and estimate B! 

#In Contrast, Generalized Least Squares:

#GLS: B = (X'C^(-1)X)^(-1)XC^(-1)y  <-- When estimating B, COVARIANCE is NOT assumed to be the identity matrix!  
#The difference between the two is the weighted component, C, of GLS and its assumptions, etc etc...

#@Brief: This section will verify all assumptions for my methodology are met before running Monte Carlo.


# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Step 3: Use Monte Carlo to Generate N Synthetic Ground Station Series
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -

#STEP by STEP: 
#       I already have C, which is set permanently. 
#       Hv = w #This is the filter used to create a vector of colored noise for a simulation.
#       w is the vector of colored noise we need for each simulation that has the same overarching noise model with 
#       varying values at each epoch.
#       H = summation(h_i) written in matrix format.
#       h_i = binom(n,k) function and it varies by station direction.
#       v is a vector with independent and identically distributed (IID) Gaussian noise
#       


#The goal of this section is to generate N Synthetic ERROR series
# TODO: all_synthetic_series is built here as a flat list via
# nested .append() calls, but fit_LS_models() below reads from it with
# all_synthetic_series[station_direction], as if it were a dict. Decide on one structure
# (e.g. dict keyed by station_direction -> list of N series) and make both ends consistent.
# Same issue applies to all_velocities / VELOCITY: this returns a list (one per direction),
# but fit_LS_models uses VELOCITY as a single scalar in the metric formula.
def generate_monte_carlo_series(H, station):
    directions = ["north", "east", "vert"]
    noise_models = []
    all_synthetic_series = []
    all_velocities = []
    for direction in directions:
        station_direction = str(station) + "_" + str(direction)
        #This sets the TRUE velocity that OLS and GLS will conform to
        VELOCITY = betas[station_direction][1].round() #I need to think about whether this is deterministic
        a, c, d, e, f = betas[station_direction][[0, 2, 3, 4, 5]]
        N = 500 #subject to change.
        size = station_lengths[station]
        all_velocities.append(VELOCITY)
        
        #Brief: This section generate N synthetic ground station motion time series.
        for simulation in range(N):
            v= np.random.normal(loc=0.0, scale=1.0, size=size)
            w = H@v 
            noise_models.append(w) #I will only begin saving persistently once I confirm the code works!
        #synthetic_series = [signal] + [noise]
        #synthetic_series = [X][B_true] + [noise]
            synthetic_series = X_matrices[station]@[a, VELOCITY, c, d, e, f] + w
            all_synthetic_series.append(synthetic_series) #I will make a dictionary instead!!

    return all_synthetic_series, all_velocities

# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Step 4: Refit OLS/GLS on each Series   
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -

#OLS ROUGH OUTLINE ONLY 
#this function will likely be in a for loop.
def fit_LS_models(station, direction, all_synthetic_series, C ,VELOCITY):
    fitted_OLS_models = {}
    station_direction = str(station) + "_" + str(direction)
    #OLS:
    # y = XB + e 
    # B = (X'X)^(-1) 
    X = X_matrices[station]
    y = all_synthetic_series[station_direction]
    OLS_betas, *_ = np.linalg.lstsq(X, y, rcond=None)
    fitted_OLS_models[station_direction] = OLS_betas

    #In Contrast, Generalized Least Squares:
    #GLS: B = (X'C^(-1)X)^(-1)XC^(-1)y 
    #TODO Very important. Cholesky Whitening. I admittedly need to revisit the math here for the following 4 lines. 
    L = np.linalg.cholesky(C)
    X_whitened = np.linalg.solve(L, X)
    y_whitened = np.linalg.solve(L, y)
    GLS_betas, *_ = np.linalg.lstsq(X_whitened, y_whitened, rcond=None)
    #@Brief: This section will recover velocity estimates for OLS and GLS from each simulation
    #AKA, create persistent storage. with open json etc etc.

    #Brief: This section will recover σ²_OLS and σ²_GLS from each simulation
    #Find residuals
    OLS_residuals = y - X@OLS_betas
    GLS_residuals = y - X@GLS_betas

    n = station_lengths[station]
    p = 6 #parameters in the model

    # TODO: double check (1) which residual vector belongs
    # in each sum of squares (both lines currently dot against the same `residuals` variable
    # from the top-level residuals.json, not GLS_residuals/OLS_residuals), and (2) whether
    # GLS's formal covariance should be built from X or from X_whitened.
    sigma_sqrd_GLS = (GLS_residuals @ residuals)/(n-p)
    Cov_GLS = sigma_sqrd_GLS * np.linalg.inv(X.T @ X)
    GLS_var_formal = Cov_GLS[1,1] 

    sigma_sqrd_OLS = (OLS_residuals @ residuals)/(n-p)
    Cov_OLS = sigma_sqrd_OLS * np.linalg.inv(X.T @ X)
    OLS_var_formal = Cov_OLS[1, 1]   # index 1 = b = velocity


    #NOTE: I am aware that I am improperly dealing with these vectors at the moment. this is a comceptual rough draft.
    GLS_cov_realism_metric = (GLS_betas[1] - VELOCITY)**2 / GLS_var_formal
    OLS_cov_realism_metric = (OLS_betas[1] - VELOCITY)**2 / OLS_var_formal #This creates z

    GLS_metric = GLS_cov_realism_metric**2
    OLS_metric = OLS_cov_realism_metric**2 #This creates z^2, which should be the Mahalanobis Distance metric from Poore
    return GLS_metric, OLS_metric
    


# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Step 5: Create Normalized σ² Vectors (Uncertainty Realism Metrics) for Comparison 
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# It is helpful to recall that z² is the Mahalanobis Distance Metric from Poore et al. at n = 1

directions = ["north", "east", "vert"]

for station in stations:
    for dir in directions:
        station_direction = str(station) + "_" + str(dir)
        C = create_parametric_C(station, dir)
        H, _, _ = create_general_power_law_cov_matrix(station, station_slopes[station_direction])
        synthetic_series, velocities_true = generate_monte_carlo_series(H, station) #TODO: Need to fix this so im not generating ALL monte carlos each time.
        GLS_metric, OLS_metric = fit_LS_models(station, dir, synthetic_series, C, velocities_true) #This pipeline is a little messy. 'll clean it up.



# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Step 6: Distribution Comparison (FOCAL POINT OF PROJECT)
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -


#@Brief: This section will plot OLS/GLS z² values on a histogram
# Compute histogram
m = _ #creating m discrete cells
GLS_hist, bin_edges = np.histogram(GLS_metric, bins=m)

print("Counts:", GLS_hist)
print("Bin Edges:", bin_edges)

plt.hist(GLS_hist, bins=5, edgecolor='black')
plt.title("Histogram of GLS Realism Metrics")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()


OLS_hist,bin_edges = np.histogram(OLS_metric, bins=m)

print("Counts:", OLS_hist)
print("Bin Edges:", bin_edges)

plt.hist(OLS_hist, bins=5, edgecolor='black')
plt.title("Histogram of OLS Realism Metrics")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()

#@Brief: This section will perform a goodness of fit test to confirm visual results
# compare the number of values in the cell compared to how many should be in the cell in a chi-sqrd dist
# compute chi sqrd test stat and compare to critical value from chi-sqr distr

# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Analysis of the significance of the OLS/GLS σ² histograms
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -




# What a doozy!

# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Future Modifications For this project
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -

# 1. Test for gaussian noise in true GNSS data
# 2. Array Procesing Extension across all stations
# 3. Utilize MLE rather than periodogram slope to define alpha 
#       (also k in some literature) to build parametric C. See Tero et al.
# 4. Pull USGS earthquake catalogue to confirm changepoints
# 5. Segment residual regimes between changepoints, then flag outliers (IQR/ Z-score)
# 6. Model diagnostics (Model-order, multicollinearity, etc)
# 7. Wrap project as real-time dynamic pipeline