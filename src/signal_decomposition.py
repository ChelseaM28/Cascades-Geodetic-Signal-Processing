# signal_decomposition.py
''' 
This script will explain and implement signal decomposition of four ground stations along the 
cascades region.
Jun 26, 2026
Version 1
Chelsea Momoh
'''


# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
#Note on project limitation: Model Structure  
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Model diagnostics are critical to increase power and reduce complexity. 
# However, in order to preserve the scope of the project, I bypass multicollinearity 
# and other predictor diagnostics.  
# With more time, I'd like to perform analyses to determine whether these are the best 
# parameters for all 4 stations. In the meantime, I make this assumption to move forward with decomposition.


# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# Import libraries and load my data from the json files
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -

import os
os.chdir("/workspaces/GNSS/data") #TODO: PRIOR TO SHOWING CODE, need to test filepath after the refactor.
import pandas as pd
import json
import numpy as np
#PSD imports
from scipy.signal import periodogram
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic


# This function was written prior tot he refactor. It is not used in current code, 
# but it is used in future code. I will soon move this directly to the covariance script.
def bin_psd(freqs, power, n_bins=30):
    # skip the zero-frequency point
    freqs = freqs[1:]
    power = power[1:]
    
    log_bins = np.logspace(np.log10(freqs.min()), np.log10(freqs.max()), n_bins)
    bin_means, bin_edges, _ = binned_statistic(freqs, power, statistic='mean', bins=log_bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    return bin_centers, bin_means

class Station:
    def __init__(self, station_id):
        print(f"Initializing station {station_id}.")
        self.directions = ["North", "East", "Vert"]
        self.station_id = station_id
        self.dataframe = pd.read_json(str(self.station_id) + ".json", orient = "records")
        self.alphas = {}
        self.betas = {}
        self.y_vectors = {}
        self.residuals = {}
        self.time_elapsed = []
        self.periodograms = {}

    def time_series_data_conversion(self):
        self.dataframe['Date'] = pd.to_datetime(self.dataframe['Date'])
        #calculating the elapsed time from epoch 0 for each row and converting to days
        #I take the date COLUMN as a series and subtract the value of the first row from EACH item in the series
        self.time_elapsed = (self.dataframe['Date'] - self.dataframe['Date'].iloc[0] ).dt.days / 365.25
        
    
    def build_design_matrix(self):
        self.X_matrix = np.column_stack(
            [np.ones(len(self.time_elapsed)),
            self.time_elapsed,
            np.sin(2*np.pi*self.time_elapsed),
            np.cos(2*np.pi*self.time_elapsed),
            np.sin(4*np.pi*self.time_elapsed),
            np.cos(4*np.pi*self.time_elapsed)
        ])
    
    def build_systems(self):
        for direction in self.directions:
            #building y vectors
            station_direction = str(direction) + "_" + str(self.station_id)
            self.y_vectors[station_direction] = self.dataframe[str(direction) + " (mm)"].values
            #building betas
            self.betas[station_direction], _, _, _ = np.linalg.lstsq(self.X_matrix, self.y_vectors[station_direction], rcond=None) 
            self.residuals[station_direction] = self.y_vectors[station_direction] - self.X_matrix @ self.betas[station_direction]

    def compute_psd_bins(self):
        for key, value in self.residuals.items():
            #I have one sample per day, and we want to define our frequency in cycles per year. 
            # let f = 365.25 (0.25 for leap years)
            freqs, power = periodogram(value, fs=365.25)
            self.periodograms[key] = (freqs, power)
            #binning
            freqs = freqs[1:]
            power = power[1:]
            log_bins = np.logspace(np.log10(freqs.min()), np.log10(freqs.max()), 30)
            bin_means, bin_edges, _ = binned_statistic(freqs, power, statistic='mean', bins=log_bins)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            mask = (bin_centers < 5) & (~np.isnan(bin_means))
            fit_freqs = bin_centers[mask]
            fit_power = bin_means[mask]
            log_f = np.log10(fit_freqs)
            log_p = np.log10(fit_power)
            slope, intercept = np.polyfit(log_f, log_p, 1)
            alpha = -slope
            self.alphas[key] = alpha 

    def plot_psd(self):
        for key, (freqs, power) in self.periodograms.items():
            plt.figure()
            plt.loglog(freqs[1:], power[1:]) #We cannot plot the (0,0 pair)
            plt.xlabel("Frequency (cycles/year)")
            plt.ylabel("Power (Variance)")
            plt.title(f"Power Spectral Density — {key}")
            plt.tight_layout()
            plt.savefig(str(key)+".png", dpi=120)
            #plt.close() keep this in mind for memory. maybe array processing or other large builds.
            

    def save_results(self):
        with open(f"{self.station_id}_alphas.json", "w") as f:
            json.dump(self.alphas, f, indent=2)
        with open(f"{self.station_id}_betas.json", "w") as f:
            json.dump({k: v.tolist() for k, v in self.betas.items()}, f, indent=2)
        with open(f"{self.station_id}_residuals.json", "w") as f:
            json.dump({k: v.tolist() for k, v in self.residuals.items()}, f, indent=2)
        with open(f"{self.station_id}_X_matrix.json", "w") as f:
            json.dump(self.X_matrix.tolist(), f, indent=2)

    def process_station(self):
        self.time_series_data_conversion()
        self.build_design_matrix()
        self.build_systems()
        self.compute_psd_bins()
        self.plot_psd()
        self.save_results()

for station_id in ["p349", "p380", "p434", "p441"]:
    station = Station(station_id)
    station.process_station()



# //// OLD EXPLANATIONS BELOW //// PRIOR TO REFACTORING ////


# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
#Step 1: Mathematical Reasoning.
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -

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

So we'ev modeled a single epoch in the east direction in the equation above, but there are thousands of lines affecting abcde&f, 
and we can't solve for that 5000+ times for each line in the dataset.
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


# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
#step 4: Characterization of noise reasoning 
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -

#@Brief: Explaining noise characterization 

'''
--Signal Decomposition--
My residual series (my residual vectors) are all time series. They are my target "signals" for analysis. 
To characterize the residuals (characterize the noise), I will decompose the signal into its 'frequency' components,
which can be understood as patterns within the data which characterize the noise.
The underlying mechanic is a discrete (due to sampling) fourier transform which 
transforms the displacement signal into frequencies. 
The fourier transform allows non-obvious patterns within the signal to be revealed before your very eyes.


Any pattern that repeats itself in this time series can be described by its frequency. 
e.g.
A pattern repeating every 10 years has 0.1 cycles per year.
A patter repeating every 6 months has 2 cycles per year. 
The lowest possible frequency for my data would have to happen about every 20 years. Or one cycle over the course of 20 years.
The Highest frequency I could detect would be HALF the sampling rate, so one cycle each 2 days.
(This is due to the effects of aliasing. At a minimum to detect a period, that period must be the length of the nyquist frequency. 


In the fourier sense, each signal is the sum of sinusoids of different frequencies. 

Power Spectral Density (PSD) is a decomposition of the signal across different frequencies (patterns).
PSD, tells, for every frequency of sinusoid, how much power the signal contains at that frequency. 

   ##Power is mathematically the average of the squared amplitude of a signal.## 

It is sometimes defined as energy per unit of time, which in this context can be analog to the 
strength of the frquency/ the strength of the noise frequency. In this case I will refer to power as variance, as 
variance is also the average of a squared quantity (avg of squared distance from the mean).
The term density comes from dividing the power spectrum by equivalent noise bandwidth (approximation of power for a given frequency) 
to lessen the effect of spectral leakage.

The plot of the PSD on a loglog plot (loglog for easier fitting) reveals the type of noise I'm dealing with.

--Types of noise--
Typically, when we estimate, for example, the velocity coefficient, the uncertainty of the 
estimate depends on the residual noise. Standard least squares assumes white (random/independent/uncorrelated) noise.
With white noise, the more data one acquires, the more accurate a prediction becomes. 
Note, I am not saying the residuals are ~N /gaussian. That is a separate assumption about the distribution of values
the residuals are drawn from. The white noise assumption is specifically about whether residuals are correlated.

However, with colored (Pink/Flicker or Random Walk) noise, noise is not independent. 
A measuremnt is affected by the last measurement. Physically, we attribute this noise to unmodeled 
physical processes. 

By plotting the PSD's shape, I can determine the type of noise I am working with, construct the 
proper covariance matrix, and plug that matrix into an appropriate uncertainty formula to get
a MORE REALISTIC ESTIMATE OF THE ERROR IN OUR VELOCITY (trend) AND ANNUAL/SEMI-ANUAL PROCESSES.

Otherwise, we'd have a pretty bad idea of how trashy our model is.

*NOTE: The model does not incorporate the comparatively large movements of the North American tectonic plate. We
are tracking the deviation of stations' movement from the plate due to seasonal loading signals, post-earthquake
deformation, etc. Recall NAM14 (the data I downloaded, see data folder) removes plate movement signals. 
This essentially 'centers data around the mean,' preventing large tectonic movement from hiding smaller signals. 
'''

#@Brief: Scipy Periodogram
'''
Periodogram takes my residual array and sampling frequency to output frequencies and power (variances)
Periodogram allows me to construct a PSD graph without the need to write out the operations done to the amplitude spectrum dervied from the fourier transform. 

Assumptions:
Periodogram requires that the mean and the variance of the signal are the same no matter where I sample. (stationary process)
For learning purposes I assume the stationarity assumption is met.
It also requires a rectangular window (that amplitude and power calculations are not affected by, for example, an attenuated signal as in a hamming window).
No attenuations or changes to the window were made, I assume the assumption is met.
It also requires a discrete spectrum. As previously mentioned, the data is descrete, not continuous. Assumption satisfied.
'''

"""
I did print(power[:10]) and got:
[2.78104284e-31 6.19735966e+01 2.85698624e+01 4.04845588e+01
 2.65258065e+00 3.20378456e+00 1.05712033e+00 4.64806175e-01
 5.13599883e-01 1.69800855e+00]
This array lists the variances of each pattern in the order corresponding to the frequencies list.
We see very high power for the first few frequencies with a large drop off near higher frequencies.
Since the power is not relatively consistent throughout the array, we know the noise is not independent/random, 
but instead colored. 

*NOTE: for the frequency array, the first number is zero, suggesting a flat-line or no pattern. A constant.
However, our model accounted for any constant term with the intercept column of X, so the periodogram will 
pick up no variance for frequencies at 0. That's why the first (freq, power) term is (0, 2.78104284e-31) (basically 0,0)

TODO: Flagging to potentially define the first number as DC. Need to look into this first.
"""

# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -# * - * - * - * - * - * 
#step 5: Characterize residual noise using Power Spectral Density (PSD) plots
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -# * - * - * - * - * - * 

# Looking at the red line, I can clearly see that there is a flattening happening at around 10^(0.7). 
# The flat section represent white noise, while the downward sloping section is colored. Maybe pink.
# NOTE: Let's set a line of demarcation at 10^(0.7), aorund 5 cycles/year

# NOTE: This shorter record for P441 has generally white noise throughout the dataset. 
#       There is not enough of a clear downward slope for polyfit to catch colored noise. Consider a Kolmogorov test without p441.


#@Brief: This section will fit the binned data for frequencies below my cutoff (for colored noise)
# For context, log(PSD)=log(A)−αlog(f). Recall we are plotting (freq, power) pairs. This function is in
# y = mx + b format, with each term having the log applied.
# α is the noise component/slope. log(f) is the input. log(A) is t  he intercept. log(PSD) is power.

#CHARACTERIZE NOISE TYPES AS FOLLOWS: white noise (α≈0), flicker (α≈1), & random walk (α≈2)
#I need to find α for each direction for each station.

#np.polyfit(x, y, 1) fits a line to the data and returns [slope, intercept]

'''
//References//
st-andrews 'spectral analysis' - https://www.st-andrews.ac.uk/~wjh/dataview/tutorials/sonogram.html
gaussian waves.com 'Power and Energy of a Signal : Demystified' - https://www.gaussianwaves.com/2013/12/power-and-energy-of-a-signal/
MATLAB 'Understanding Power Spectral Density and the Power Spectrum' - https://www.youtube.com/watch?v=pfjiwxhqd1M
'''