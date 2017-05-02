import numpy as np
import matplotlib.pyplot as plt
import random



npts = 100000
mu, sigma = 10, 1
x = mu - sigma*np.random.chisquare(mu, size=npts)

mean   = np.mean(x)
median = np.median(x)
rms    = np.sqrt(np.mean(np.square(x)))

plt.figure()

lw = 2
plt.hist(x, bins=50, alpha=0.4)
plt.axvline(mean, color='r', label='Mean', lw=lw)
plt.axvline(median, color='g', label='Median', lw=lw)
plt.axvline(rms, color='b', label='RMS', lw=lw)
plt.legend()
plt.show()
