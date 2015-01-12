import numpy as np
import pylab

mu, sigma = 0.0076, 0.0076*0.1/2
Wt = mu + sigma * np.random.randn(2000)
pylab.hist(Wt,20)
pylab.show() 
