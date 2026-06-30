# signal_decomposition.py
''' 
This script will explain and implement signal decomposition of four ground stations along the 
cascades region.
Jun 26, 2026
Version 1
Chelsea Momoh
'''

#Step 1: Mathematical Reasoning.

'''
Our data consists of three vectors I will name y_north, y_east, and y_vertical. 
They are the displacement in the north, east, and vertical directions at each epoch.

I found/assume each displacement, y, can be approximated using 
    - tectonic velocity
    - annual amplitude and phase (wiggle caused by snow melt/groundwater)
    - semiannual amplitude and phase (wiggle caused by other physcial processes)

For example, for a single epoch, we represent a displacement in the east direction as such:

y_east_i = a + bt_i + csin(2PIt_i) + dcos(2PIt_i) + esin(4PIt_i) + fcos(4PIt_i) , 

Where 
a -> value of intercept @ t = 0
b -> linear coefficient of tectonic velocity (tectonic velocuty is constant)
c & d -> used in combination as linear coefficients to weight the amplitude and phase of the annual cycle
e & f -> used in combination as linear coefficients to weight the amplitude and phase of the SEMI-annual cycle
t -> decimal years since reference point (this is the time)
y -> millimiters of displacement in the east direction at the ith epoch


We don't represent the quarterly cycle to avoid overfitting. Also because quarterly cycles don't have significant physical processes.

So we'ev modeled a single epoch in the east direction in the equation above, but there are thousands of lines affecting abcde&f, and we can't solve for that 5000+ times for each line in the dataset.
So instead, we solve all at once using a matrix formula.

XB_east = Y_east

where 
X = [1 t_1 sin(2PIt_1) cos(2PIt_1) sin(4PIt_1) cos(4PIt_1)
     1 t_2 sin(2PIt_2) cos(2PIt_2) sin(4PIt_2) cos(4PIt_2)
     1 t_3 sin(2PIt_3) cos(2PIt_3) sin(4PIt_3) cos(4pPIt_3)
     . . .
     .
     .
     1 t_n sin(2PIt_n) cos(2PIt_n) sin(4PIt_n) cos(4PIt_n) ]

B = transpose([a b c d e f]) for the east direction

So X is the same for all three directions. 
But Y_north, Y_east, and Y_vertical change, so we'll need different coefficients. In other words, 
different B_east, B_vertical, B_north.

Got me?
'''


#Step 2: Load my data from the json files
import os
os.chdir("/workspaces/GNSS/data")
import pandas as pd
import json
import numpy as np

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



#step 3: Fit the design models

#@Brief: This section will utilize datetime to perform datetime arithmetic and find t. 
#The specific dataset I use is not arbitrary. Each station dataset has a different number of elapsed days, so i have to repeat this process four times.
# The formula to convert the number of days to decimal years and get time deltas (changes in time) is given as:

elapsed = p349['Date'] - p349['Date'].iloc[0] 
t_p349 = elapsed.dt.days / 365.25

elapsed = p380['Date'] - p380['Date'].iloc[0] 
t_p380 = elapsed.dt.days / 365.25

elapsed = p434['Date'] - p434['Date'].iloc[0] 
t_p434 = elapsed.dt.days / 365.25

elapsed = p441['Date'] - p441['Date'].iloc[0] 
t_p441 = elapsed.dt.days / 365.25


#@Brief: This section will construct the design matrices for X for all four stations.

X_p349 = np.column_stack([
    np.ones(len(t_p349)),  # ones
    t_p349,  # t
    np.sin(2 * np.pi * t_p349),  # sin(2π t)
    np.cos(2 * np.pi * t_p349),  # cos(2π t)
    np.sin(4 * np.pi * t_p349),  # sin(4π t)
    np.cos(4 * np.pi * t_p349)   # cos(4π t)
])

X_p380 = np.column_stack([
    np.ones(len(t_p380)),  # ones
    t_p380,  # t
    #I must use np.sin for element-wise operations. Numpy is great for element-wise stuffs!
    np.sin(2 * np.pi * t_p380),  # sin(2π t)
    np.cos(2 * np.pi * t_p380),  # cos(2π t)
    np.sin(4 * np.pi * t_p380),  # sin(4π t)
    np.cos(4 * np.pi * t_p380)   # cos(4π t)
])

X_p434 = np.column_stack([
    np.ones(len(t_p434)),  # ones
    t_p434,  # t
    np.sin(2 * np.pi * t_p434),  # sin(2π t)
    np.cos(2 * np.pi * t_p434),  # cos(2π t)
    np.sin(4 * np.pi * t_p434),  # sin(4π t)
    np.cos(4 * np.pi * t_p434)   # cos(4π t)
])

X_p441 = np.column_stack([
    np.ones(len(t_p441)),  # ones
    t_p441,  # t
    np.sin(2 * np.pi * t_p441),  # sin(2π t)
    np.cos(2 * np.pi * t_p441),  # cos(2π t)
    np.sin(4 * np.pi * t_p441),  # sin(4π t)
    np.cos(4 * np.pi * t_p441)   # cos(4π t)
])

#@Brief: This section will extract the y vectors from the dataframes and perform 
# least squares approx to solve for B for all four stations.

#This gives me the y vectors for each direction for each station 
y_north_p349 = p349['North (mm)'].values
y_east_p349 = p349['East (mm)'].values
y_vert_p349 = p349['Vertical (mm)'].values


y_north_p380 = p380['North (mm)'].values
y_east_p380 = p380['East (mm)'].values
y_vert_p380 = p380['Vertical (mm)'].values

y_north_p434 = p434['North (mm)'].values
y_east_p434 = p434['East (mm)'].values
y_vert_p434 = p434['Vertical (mm)'].values

y_north_p441 = p441['North (mm)'].values
y_east_p441 = p441['East (mm)'].values
y_vert_p441 = p441['Vertical (mm)'].values


#@Brief: this next section gives me the coefficients for each direction for all four stations.
#The underscores get rid of unneeeded data (e.g., the residuals provided here are not a time series, but a summary stat, not helpful.)
beta_north_p349, _, _, _ = np.linalg.lstsq(X_p349, y_north_p349, rcond=None) 
beta_east_p349, _, _, _ = np.linalg.lstsq(X_p349, y_east_p349, rcond=None) 
beta_vert_p349, _, _, _ = np.linalg.lstsq(X_p349, y_vert_p349, rcond=None) 

beta_north_p380, _, _, _ = np.linalg.lstsq(X_p380, y_north_p380, rcond=None) 
beta_east_p380, _, _, _ = np.linalg.lstsq(X_p380, y_east_p380, rcond=None) 
beta_vert_p380, _, _, _ = np.linalg.lstsq(X_p380, y_vert_p380, rcond=None) 

beta_north_p434, _, _, _ = np.linalg.lstsq(X_p434, y_north_p434, rcond=None) 
beta_east_p434, _, _, _ = np.linalg.lstsq(X_p434, y_east_p434, rcond=None) 
beta_vert_p434, _, _, _ = np.linalg.lstsq(X_p434, y_vert_p434, rcond=None) 

beta_north_p441, _, _, _ = np.linalg.lstsq(X_p441, y_north_p441, rcond=None) 
beta_east_p441, _, _, _ = np.linalg.lstsq(X_p441, y_east_p441, rcond=None) 
beta_vert_p441, _, _, _ = np.linalg.lstsq(X_p441, y_vert_p441, rcond=None) 

'''
I wrote print(beta_north_p349) and got
[-3.33226909  6.91818185  0.03229881 -0.7911755   0.26515507 -0.06557577]
Notice a = -3.33 despite the data showing a displacement of zero at this point. This is just noise around the intercept!
'''

print("Completed building design matrices.")

#@Brief: This section will find the residuals for each groundstation's directions.
#Residuals are milimeters of variation the model did not explain.
#residuals=y−Xβ

resid_p349_north = y_north_p349 - X_p349@beta_north_p349
resid_p349_east = y_east_p349 - X_p349@beta_east_p349
resid_p349_vert = y_vert_p349 - X_p349@beta_vert_p349

resid_p380_north = y_north_p380 - X_p380@beta_north_p380
resid_p380_east = y_east_p380 - X_p380@beta_east_p380
resid_p380_vert = y_vert_p380 - X_p380@beta_vert_p380

resid_p434_north = y_north_p434 - X_p434@beta_north_p434
resid_p434_east = y_east_p434 - X_p434@beta_east_p434
resid_p434_vert = y_vert_p434 - X_p434@beta_vert_p434

resid_p441_north = y_north_p441 - X_p441@beta_north_p441
resid_p441_east = y_east_p441 - X_p441@beta_east_p441
resid_p441_vert = y_vert_p441 - X_p441@beta_vert_p441

print(f"Finished building residual series.")

#step 4: Characterization of noise reasoning 

'''
My residual series (my residual vectors) are all time series. 
(In the Fourier sense, they are the sum of sinusoids.)
Any pattern that repeats itself in this time series can be described by its frequency. 
e.g.
A pattern repeating every 10 years has 0.1 cycles per year.
A patter repeating every 6 months has 2 cycles per year. 
The lowest possible frequency for my data would have to happen every 20 years. Or one cycle over the course of 20 years.
The Highest frequency would equal the sampling rate, so it would happen each day. One cycle each 2 days.

Power Spectral Density (PSD) tells, for every frequency of sinusoid, how much power (variance) the signal
contains at that frequency. PSD is a decomposition of the signal across different frequencies (patterns).


The plot of the PSD on a loglog plot (loglog for easier fitting) reveals the type of noise I'm dealing with.

Typically, when we estimate, for example, the velocity coefficient, the uncertainty of the 
estimate depends on theresidual noise. Standard least squares assumes white (random/normal/gaussian) noise.
With white noise, the more data one acquires, the more accurate a prediction becomes. 

However, with colored (Pink/Flisker or Random Walk) noise, noise is not independent. 
Each measuremnt is affected by the last. Physically, we attribute this noise to unmodeled 
physical processes.

By plotting the PSD's shape, I can determine the type of noise I am working with, construct the 
proper covariance matrix, and plug that matrix into an appropriate uncertainty formula to get
a MORE REALISTIC ESTIMATE OF THE ERROR IN OUR VELOCITY (trend) AND ANNUAL/SEMI-ANUAL PROCESSES.

Otherwise, we'd have a pretty bad idea of how trashy our measurements are.

*NOTE: We are not tracking the large signals of the North American tectonic plate's movement. We
are tracking the deviation of stations' movement from the plate due to seasonal loading signals, post-earthquake
deformation, etc. Recall NAM14 (the data I downloaded, see data folder) removes plate movement signals. 
This essentially 'centers data around the mean,' preventing large tectonic movement from hiding smaller signals. 

Ya feel?  

'''


#step 5: Characterize residual noise using Power Spectral Density (PSD) plots






