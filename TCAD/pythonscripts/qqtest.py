from scipy import stats
import matplotlib.pyplot as plt
import numpy as np

nsample = 100


ax1 = plt.plot(1)
x = stats.t.rvs(3, size=nsample)
res = stats.probplot(x, plot=plt)




plt.show()
