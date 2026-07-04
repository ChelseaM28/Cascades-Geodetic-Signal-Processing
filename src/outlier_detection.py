# outlier_detection.py
''' 
This script will implement IQR, Z-score flagging, and a change-point detection model (PELT).
Jun 30, 2026
Version 1
Chelsea Momoh
'''

#Step 1: Load libraries and data
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
#Step 2: Outlier Detection Reasoning

'''
My residuals, the difference between my model and actual data, shoud be explained as noise, whether white or colored. 
However, outliers will occur under a set of circumstances, and these events should be flagged:
- Equipment changes/malfunctions, earthquakes, atmospheric noise, and glitches could afffect daa quality and lead to an outlier.
This script will utilize two tools, IQR/Z-score flaggin and change-point detection (PELT) to flag outliers.

Source: https://pipiras.sites.oasis.unc.edu/
Outlier: "have only an instantaneous effect" on data quality
Change points: "effect decays over time.... sustained for the entire series (or an extended portion of the series)"
//end source

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

    Source: https://www.lancaster.ac.uk/~romano/teaching/2425MATH337/4_algos_and_penalties.html
    PELT (Pruned Exact Linear Time): In the past, to segment data according to change points, binary segmentation and optimal 
partitioning were used. However, BS is prone to error and OP is computationally expensive. PELT solves these issues by
"reduc[ing] the [number] of [changepoint] checks to be performed at each iteration," also called pruning. It introduces a penalty
to discourage the algorithm from continuosly adding changepoints (overfitting).
PELT
INPUT: Time series, penalty 
OUTPUT: Optimal changepoint vector 


'''


#Step 3: Implement PELT Change point model first so I can segment data for IQR and Z-score flagging

#@Brrief: This section will test the PELT algorithm outputs.
'''signal = residuals['p349_north']  

algo = rpt.Pelt(model="rbf").fit(signal)
breaks = algo.predict(pen = 10)

I wrote:
print("Change points detected at indices:", breaks)
print(len(residuals['p349_north']))
and got:
7524
Change points detected at indices: [515, 895, 1565, 2695, 3485, 3770, 4415, 5235, 5600, 6125, 6375, 6855, 7070, 7524]
These are indices in the residual time series, in other words, these are EPOCHS at which a change point occurs.
So day 515, day 895, day 1565, etc. NOT frequencies.
Note, the last value is just the end of the time series, not a change point.
'''

#@Brief: In this section, I will look closely at the penalty values and determine whether it is reasonable.
'''
changepoints = {}
for key, value in residuals.items():
    signal = value
    algo = rpt.Pelt(model="rbf").fit(signal)
    breaks = algo.predict(pen=10)
    
    station = key.split("_")[0]
    dates = station_dates[station]
    changepoints[key] = [dates[i] for i in breaks[:-1]]  # Exclude the last index which is the end of the series
    print(f"Change points detected for {key} at dates:", changepoints[key])

My initial penalty is way too permissive. 
Source: https://academic.oup.com/gji/article/204/1/480/635055?login=false
With at least 700 days (but usually more) of data, a given station might have 1.8 changepoints on average. 
49% caused by documented equipment changes, 31% by earthquakes, and 20% due to unknown causes." --> [should discuss anomaly detection in context of space domain awareness in final write-up]
NOTE: Next, I will experiment with penalty values and create charts and persistent storage of results.
'''

#The length of the dataset causes my computer to crash. 
#To fix this, I run the PELT algo using a higher jump number (from baseline 5 to 10), meaning it will have fewer checks for changepoints. 
#I'm beginning to see why computer scientists are all about efficiency/bestcase/worstcase scenarios.

#@Brief: In this section, I determine the best range of penalty values (for a single station) before constructing a loop to run PELT on the residuals of all stations. 
#Penalty values are unique to each station, but it is helpful to have a range of values to test.
'''print(f"Length of residuals for p349_north: {len(residuals['p349_north'])}")
signal = residuals['p349_north']
pen_values = [110,112,114,116,118,120]
counts = []
algo = rpt.Pelt(model="rbf", jump=10).fit(signal)

for pen in pen_values:
    print("Beginning penalty tracking: ", pen)
    breaks = algo.predict(pen=pen)
    counts.append(len(breaks) - 1)

plt.figure()
#will try linear first, not loglog
plt.plot(pen_values, counts, color='red', linewidth=2, marker='o', label="penalty vs. # of change points")
plt.xlabel("Penalty Value")
plt.ylabel("Change Points Detected")
plt.title("PenaltiesXChange Points p349_north")
plt.suptitle(f"jump={10}, model=rbf", fontsize=9, y=0.93)
plt.legend()
plt.tight_layout()
plt.savefig("p349_north_penalties.png", dpi=120)
print("Saved p349_north_penalties")
'''


#@Brief: This section will loop through all stations and create penalty plots
pen_values = [110,112,114,116,118,120]

for key, value in residuals.items(): #I always forget to add .items() !!
    signal = value
    jump = 20
    algo = rpt.Pelt(model="rbf", jump=jump).fit(signal) #Adding larger jumps because I don't want the runtime to be too long.
    counts = []
    print(f"Drawing penalties for {key}")
    for pen in pen_values:
        print("Beginning penalty tracking: ", pen)
        breaks = algo.predict(pen=pen)
        counts.append(len(breaks) - 1)
    plt.figure()
    plt.plot(pen_values, counts, color='red', linewidth=2, marker='o', label = "Penalty vs. Detected Change Points")
    plt.xlabel("Penalty Value")
    plt.ylabel("Change Points Detected")
    plt.title("Penalties X Change Points for " + key)
    plt.suptitle(f"jump={jump}, model=rbf", fontsize=9, y=0.93)
    plt.legend()
    plt.tight_layout()
    plt.savefig(key + "_penalties.png", dpi=120)
    print("Saved " + key + "_penalties.png")
print("Finished drawing penalty plots for all stations")

#I will fill out this dictionary by eye-balling the penalty plots
final_penalties = {"p349_north": 119, "p349_east": 115, "p349_vert": None,
                   "p380_north": None, "p380_east": None, "p380_vert": None,
                   "p434_north": None, "p434_east": None, "p434_vert": None,
                   "p441_north": None, "p441_east": None, "p441_vert": None,
                   } #NOTE: Most graphs were flat lines. I need to implement a broader search. Left off here.

#Step 4: Implement IQR and Modified Z-score Flagging  