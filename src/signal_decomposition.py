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

Where a -> value of intercept @ t = 0
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


#this next section gives me the coefficients for each direction for all four stations.
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

print("Completed building. design matrices.")


#step 4: Characterize residuals
