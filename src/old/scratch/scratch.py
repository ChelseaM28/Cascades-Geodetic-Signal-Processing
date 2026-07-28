#This file stores dead code

#NOTE From changepoint_detection.py

'''
This was me testing the lower triangular matrix in a collab notebook
import numpy as np
h = [1]
k = -2 #doesn't work unless k is negative!
n = 10
for i in range(n):
  index = 1 + i
  h.append((index - (k/2) - 1)*(h[index-1]/index))
print(h)
H = np.zeros((n,n)) #I am going to convert h into matrix form, H
for i in range(n): 
      #h[h::-1]  = h[start:stop:skip]
      H[i, :i+1] = h[i::-1] #The whole list gets filled in!
print(H)

'''


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


#NOTE From signal decomposition

'''
#I'm only commenting this out so I don't continue to create more plots each time I run!
for key, (freqs, power) in PSD_set.items():
    plt.figure()
    plt.loglog(freqs[1:], power[1:]) #We cannot plot the (0,0 pair)
    plt.xlabel("Frequency (cycles/year)")
    plt.ylabel("Power")
    plt.title(f"PSD — {key}")
    plt.tight_layout()
    plt.savefig(str(key)+".png", dpi=120)'''


'''
#Commenting this out to avoid further images!
freqs, power = periodogram(residuals['p349_north'], fs=365.25)
bin_centers, bin_means = bin_psd(freqs, power)

plt.figure()
plt.loglog(freqs[1:], power[1:], alpha=0.3, label="raw")
plt.loglog(bin_centers, bin_means, color='red', linewidth=2, label="binned mean")
plt.xlabel("Frequency (cycles/year)")
plt.ylabel("Power")
plt.title("PSD — p349_north (binned)")
plt.legend()
plt.tight_layout()
plt.savefig("p349_north_binned.png", dpi=120)'''


'''
I wrote print(beta_north_p349) and got
[-3.33226909  6.91818185  0.03229881 -0.7911755   0.26515507 -0.06557577]
Notice a = -3.33 despite the data showing a displacement of zero at this point. This is just noise around the intercept!
'''


'''
#commenting out a plot again
freqs, power = periodogram(residuals['p441_east'], fs=365.25)
bin_centers, bin_means = bin_psd(freqs, power)

plt.figure()
plt.loglog(freqs[1:], power[1:], alpha=0.3, label="raw")
plt.loglog(bin_centers, bin_means, color='red', linewidth=2, marker='o', label="binned mean")
plt.axvline(5, color='gray', linestyle='--', label="cutoff (5 cyc/yr)")
plt.xlabel("Frequency (cycles/year)")
plt.ylabel("Power")
plt.title("PSD — p441_east (binned)")
plt.legend()
plt.tight_layout()
plt.savefig("p441_east_binned.png", dpi=120)'''


'''
Finishing off the rest of the plots.
already_plotted = {"p349_north", "p441_east"}
 
for key, value in residuals.items():
    if key in already_plotted:
        continue
 
    freqs, power = periodogram(value, fs=365.25)
    bin_centers, bin_means = bin_psd(freqs, power)
 
    plt.figure()
    plt.loglog(freqs[1:], power[1:], alpha=0.3, label="raw")
    plt.loglog(bin_centers, bin_means, color='red', linewidth=2, marker='o', label="binned mean")
    plt.axvline(5, color='gray', linestyle='--', label="cutoff (5 cyc/yr)")
    plt.xlabel("Frequency (cycles/year)")
    plt.ylabel("Power")
    plt.title(f"PSD — {key} (binned)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{key}_binned.png", dpi=120)
    plt.close()
 
print("Done plotting remaining PSDs.")'''

