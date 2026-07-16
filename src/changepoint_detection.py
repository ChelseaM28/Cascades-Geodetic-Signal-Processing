# changepoint_detection.py
''' 
This script will implement a change-point detection model (PELT). 
Jun 30, 2026
Version 1
Chelsea Momoh
'''

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
# Step 2: Outlier Detection Reasoning
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
'''
My residuals, the difference between my model and actual data, shoud be explained as noise, whether white or colored. 
However, outliers will occur under a set of circumstances, and these events should be flagged:
- Equipment changes/malfunctions, earthquakes, atmospheric noise, and glitches could afffect data quality and lead to an outlier.
This script will utilize two tools, IQR/Z-score flagging and change-point detection (PELT) to flag outliers.

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

    IQR would flag values outside a given range: [Q1 - 1.5*IQR, Q3 + 1.5*IQR] and does not assume Gaussian noise.
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

# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - *
# Note on project limitation: Implementing IQR and Modified Z-score Flagging  
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - *
# Whether using OLS or GLS, outliers negatively impact model inferences and uncertainty. Identifying 
# outliers and refitting data is a key part of statistical analysis. However, to preserve the scope 
# of the project, I bypass full data refitting. The effects of this mean GLS uncertainty formula is
# validated with the assumed noise model, not a corrected noise model. With more time, I would certainly 
# include additional outlier detection. 



# * - * - * - * - * - * - * - * - * - * - * 
# Step 3: Implement PELT Change point model 
# * - * - * - * - * - * - * - * - * - * - * 

#@Brief: This section will test the PELT algorithm outputs.
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
def pen_testing():
    pen_values = []

    for key, value in residuals.items(): #I always forget to add .items() !!
        if key == 'p434_east' or key == 'p434_vert' or key =='p441_east':
            pen_values = np.arange(14, 80, 15)
            jump = 20
        else:
            continue
        
        signal = value
        
        #model = rbf is 'radial basis function,' a model used to detect changepoints without white noise assumptions.
        #The cost function used to detect changepoints is based on the similarity between groups of points rather than
        #the squared difference between terms that least squares often utilizes.
        #rbf (k(x_i, x_j) = exp( −||x_i − x_j||² / σ² )) will output values between 0 and 1 and compare similarity.

        #One issue with colored noise is changepoint detection might recognize walk noise as slow long term drift and 
        #classify noise as a changeopint erroneously. rbf, as useful as it is, doesn't distinguish between drift caused by
        #a physical process (like an earthquake) and the causeless drift of flicker/random walk. It seems to take human judgement 
        #to differentiate between the two. Maybe that is why automation of change points is so difficult.
        #I'd like to explore any resources giving ground truth about genuine changepoints. 

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
        plt.savefig(key + "round_4_penalties.png", dpi=120)
        print("round_4_Saved " + key + "_penalties.png")
    print("Finished drawing penalty plots for all stations")

#I will fill out this dictionary using visual representations: the penalty plots
final_penalties = {"p349_north": 119, "p349_east": 115, "p349_vert": 140,
                   "p380_north": 115, "p380_east": 115, "p380_vert": 55,
                   "p434_north": 65, "p434_east": 48, "p434_vert": 44,
                   "p441_north": 185, "p441_east": 29, "p441_vert": 14,
                   } 


'''
I use the elbow method to select penalty scores for PELT.
Elbow method: I look for a steep drop or a plateau surrounded by steep drops. I choose a penalty score 
in the drop or within the plateu. This is to prevent permissive penalty scores (high number of detected changepoints suggest overfitting)
and also prevents overly strict penalty scores (possibly missing out on real changepoints)

The following are the ideal penalty ranges for each station based on penalty plots. 

#--#ROUND 2#--#
p349_vert: between 130 - 150 penalty score (stable plateu)

p380_north: between 100 - 125 (stable plateu)
p380_east: between 105 - 125 (steep drop)
p380_vert: wide sweep again with range less less than 100. Ideally 50 - 100 (steep drop at beginning of graph) 

p434_north: wide sweep again with range less less than 100. Ideally 50 - 100 (steep drop at beginning of graph) 
p434_east: wide sweep again. Completely FLAT LINE. 0 changepoints detected
p434_vert: wide sweep again. Completely FLAT LINE. 0 changepoints detected

p441_north: between 180 - 200 (steep drop)
p441_east: wide sweep again. Completely FLAT LINE. 0 changepoints detected
p441_vert: wide sweep again. Completely FLAT LINE. 0 changepoints detected


#--#ROUND 3#--#

p380_east: 115 (flatline @ 2 changepoints)
p380_vert: 55 (flatline @ 3 changepoints followed by steep drop from earlier)

p434_north: 65 (flatline 1t 2 changepoints from 50 to 80)
p434_east: wide sweep again at np.arange(15, 60, 10)
p434_vert: wide sweep again at np.arange(14, 80, 15)

p441_north: 185
p441_east: wide sweep again at np.arange(14, 80, 15)
p441_vert: 14

#--#ROUND 4#--#
p434_east: 48
p434_vert: 44
p441_east: 35

'''

#@Brief: This section will run PELT changepoint detection and align changepoints 
#        with dates in each station's timeseries.

changepoints = {}
def pelt():
    print("-- Beginning Changepoint Detection --")
    for key, value in residuals.items():
        signal = value
        algo = rpt.Pelt(model="rbf", jump = 20).fit(signal)
        breaks = algo.predict(pen=final_penalties[key])
        
        station = key.split("_")[0]
        dates = station_dates[station]
        changepoints[key] = [str(dates[i]) for i in breaks[:-1]]  # the date format wont work with json, but str will
        print(f"Change points detected for {key} at dates:", changepoints[key])
    
    with open("changepoints.json", "w") as f:
        json.dump(changepoints, f, indent=4)
    print("-- Completed json changepoint records! --")

def run_single_station(key, pen, jump=20):
    signal = residuals[key]
    algo = rpt.Pelt(model="rbf", jump=jump).fit(signal)
    breaks = algo.predict(pen=pen)

    station = key.split("_")[0]
    dates = station_dates[station]
    changepoints[key] = [str(dates[i]) for i in breaks[:-1]]
    print(f"Change points detected for {key} at dates:", changepoints[key])

    with open("changepoints.json", "w") as f:
        json.dump(changepoints, f, indent=4)
    print("-- Updated changepoints.json --")

with open("changepoints.json", "r") as f:
    changepoints = json.load(f)

#Previously I chose the wrong penalty score- 35, which was too strict. 
# i need to re-run with pen = 29, according to my penalty plot from round 4
run_single_station('p441_east', pen=29)

# NOTE: I set the jump value to a permanent 20. a small jump value was necessary 
# to detect certain changepoints. where the penalty needed to be lower. But, 
# realistically, computationally, having a jump size of 5 or 10 isn't realistic. 
# I know the penalty score that is resonable for each station, which is arguably 
# more important because the jump affects how often the model conducts the search, 
# not whether the model will throw out a potential candidate due to an overly 
# strict or permissive penalty. 
# Actually, with more computational resources and time, I would like to have 
# considered how long it takes the noise resulting from an earthquake to decay to 
# ensure my jump values aren't too high. PELT has jump-limited imprecision, but my 
# later per-epoch outlier detection methods will detect outliers within each regime
# that the PELT may not have been able to catch.


