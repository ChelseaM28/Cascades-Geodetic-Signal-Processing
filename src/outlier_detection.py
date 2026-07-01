# outlier_detection.py
''' 
This script will implement IQR, Z-score flagging, and a change-point detection model (PELT).
Jun 30, 2026
Version 1
Chelsea Momoh
'''

#Step 1: Load libraries and data
import os
os.chdir("/workspaces/GNSS/data")
import pandas as pd
import numpy as np  
import json

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


#Step 2: Outlier Detection Reasoning

'''
My residuals, the difference between my model and actual data, shoud be explained as noise, whether white or colored. 
However, outliers will occur under a set of circumstances, and these events should be flagged:
- Equipment changes/malfunctions, earthquakes, atmospheric noise, and glitches could afffect daa quality and lead to an outlier.
This script will utilize two tools, IQR/Z-score flaggin and change-point detection (PELT) to flag outliers.

Source: https://pipiras.sites.oasis.unc.edu/
Outlier: "have only an instantaneous effect" on data quality
Change points: "effect decays over time.... sustained for the entire series (or an extended portion of the series)"

Causes of outliers include atmospheric interruption or known noisy GNSS signals
Causes of change points include equipment changes (antenna or receiver swapped out) or earthquakes
--> The velocity of the signal after an earthquake changes for a period of time after the earthquake until it decays to baseline 

IQR/Z-score flagging catches outliers while PELT detection catches change points.

Methodology:
    Z-Score measure how many standard deviations a data point is fromt he mean. However it assumes Gaussian noise.
In our case, we are measuring how far the residual is from the mean.
z = (x - mean) / std

    IQR will flag values outside a given range: [Q1 - 1.5*IQR, Q3 + 1.5*IQR] and does not assume Gaussian noise.
However, it is less sensitive to outliers than Z-score (massive outliers will affect standard deviation more than a quartile).

    PELT (Pruned Exact Lenear Time) NOTE: CONTINUE HERE TOMORROW



'''


#Step 3: Implement IQR and Modified Z-score Flagging 


#Step 4: Implement PELT Change point model