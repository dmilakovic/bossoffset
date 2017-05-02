import numpy as np
import matplotlib.pyplot as p
import pandas as pd

pd.set_option('display.max_rows', 500)
pd.set_option('display.width', 120)

fig, axs = p.subplots()

quantiles = [0.025, 0.16, 0.5, 0.84, 0.975]

# generate random data
ncols = 10
nrows = 100
df1 = pd.DataFrame(np.random.randn(nrows,ncols))

# randomly choose positions to fill with nans
nn = 150
for i in xrange(nn):
    c = np.random.randint(low = 0, high = ncols, size = 1)
    r = np.random.randint(low = 0, high = nrows, size = 1)
    df1.loc[r,c] = np.nan

# calculate medians using pandas
pandasq = df1.quantile(q=0.5,axis='columns').transpose()

# calculate medians using numpy
numpyq = np.zeros(shape=(nrows,))
for i in xrange(nrows):
    numpyq[i] = np.nanpercentile(df1.iloc[i], 50)
numpyq = pd.Series(numpyq)

# create dataframe with medians
df2 = pd.DataFrame({'pandas':pandasq, 'numpy':numpyq})

# plot medians
df2.plot(ax=axs)
axs.set_ylim(-1,1)
p.legend()
p.show()
