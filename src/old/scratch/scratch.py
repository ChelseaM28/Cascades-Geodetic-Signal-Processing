#This file stores dead code

#From signal decomposition

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