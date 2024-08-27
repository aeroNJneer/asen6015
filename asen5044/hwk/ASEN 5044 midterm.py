#ASEN 5044 - HWK #4

# 2.16 Generate 10,000 samples of (x1 + x2)/2, where each 2% is a random number 
# uniformly distributed on [-1/2, +1/2]. Plot the 50-bin histogram. Repeat for 
# (x1+ x2 + x3 + x4)/4. Describe the difference between the two histograms


import numpy as np
import random
import matplotlib.pyplot as plt
import pandas as pd



n=10000
rand_percent = 0.02
rangelow = -1/2
rangehigh = 1/2
bin = 50

#fxn1 = (x1 + x2)/2
#fxn2 = (x1 + x2 +x3 + x4)/4
# 
x1 = []
x2 = []
fxn2 = []

for i in range(0,1000):
    x1.append(random.random())
    x2.append(random.random())
    fxn2.append((x1[i] + x2[i])/2)
    
    #float x4 = random.sample(range(rangelow,rangehigh),10000)
#float x4 =  random.sample(range(rangelow,rangehigh),10000)


plt.hist(fxn2,bins=50)



