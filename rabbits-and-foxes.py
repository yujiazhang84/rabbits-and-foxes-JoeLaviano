# In[1]:

# Kinetic Monte Carlo Method

import numpy as np
from matplotlib import pyplot as plt
import random

# #import cProfile
# import time
# import pstats

# Intializing Constants
k1 = 0.015  # days^-1 Growth of Rabbits
k2 = 0.00004  # days^-1foxes^-1 Rabbits eaten by foxes
k3 = 0.0004  # days^-1 foxes eating rabbits
k4 = 0.04  # days^-1 Death of foxes
R0 = 400  # Initial Rabbits
F0 = 200  # Initial Foxes
days = 600

# Initialize Counter for Extinction Case
FDeathCount = 0
RDeathCount = 0
trials = 1000

# Initialize Arrays for Second Peak Calculations
SecondPeak = []
TimeSecondPeak = []

for i in range(trials):

    # Reset for Each Trial
    rabbits = R0
    foxes = F0
    t = [0]
    KMCR = [rabbits]
    KMCF = [foxes]

    while t[-1] < days:

        # Rates multiply by the last item in the array (Most Recent Population)
        Rb = k1 * rabbits
        Rd = k2 * rabbits * foxes
        Fb = k3 * rabbits * foxes
        Fd = k4 * foxes
        TotalRate = Rb + Rd + Fb + Fd

        # Generate Random Number
        Rando = np.random.rand() * TotalRate

        if rabbits == 0:
            # Rabbits are Dead
            RDeathCount = RDeathCount + 1
            FDeathCount = FDeathCount + 1
            break
        elif foxes == 0:
            # Foxes are Dead
            FDeathCount = FDeathCount + 1
            break
        if Rando < Rb:
            # Rabbit is Born
            rabbits = rabbits + 1
            KMCR.append(rabbits)
            KMCF.append(foxes)

        elif Rando < (Rb + Rd):
            # Rabbit Dies
            rabbits = rabbits -1
            KMCR.append(rabbits)
        elif Rando < (Rb + Rd + Fb):
            # Fox is Born
            foxes = foxes + 1
            KMCR.append(rabbits)
            KMCF.append(foxes)
        elif Rando < TotalRate:
            # FoxDies
            foxes = foxes - 1
            KMCR.append(rabbits)
            KMCF.append(foxes)
        #use random.expovariate to speed up and move it out of the loop
        t.append(t[-1] + random.expovariate(TotalRate))

    # Getting Second Peaks and Second Peak Times for each Trial
    KMCSecondHalf = KMCF[len(KMCF) // 2:]
    TSecondHalf = t[len(KMCF) // 2:]
    Peak2 = np.max(KMCSecondHalf)
    SecondPeak.append(Peak2)
    Timeslot = KMCSecondHalf.index(Peak2)
    TimeSecondPeak.append(TSecondHalf[Timeslot])

# I tried converting list to array here, but it didn't help with speeding up...
#    arrayF = np.array(KMCF)
#    arrayT = np.array(t)
#    HalfArrayF = np.array_split(arrayF, 2)
#    HalfArrayT = np.array_split(arrayF, 2)
#    KMCFSecondHalf = HalfArrayF[1]
#    TSecondHalf = HalfArrayT[1]
#    Peak2 = np.max(KMCFSecondHalf)
#    Peak2Location = TSecondHalf[KMCFSecondHalf.argmax()]
#    SecondPeak.append(Peak2)
#    TimeSecondPeak.append(Peak2Location)

    # Stole plotting on the same graph from @Justin
#    plt.plot(t,KMCR,'b')
#    plt.plot(t,KMCF,'r')
#    plt.legend(loc='best')

# plt.xlabel('Days')
# plt.ylabel('Population')
# plt.show()
print("Foxes died out ", FDeathCount, " times and Rabbits died out ", RDeathCount, " times.")

# In[2]:

# AvgSecondPeak = np.average(SecondPeak)

# My way of getting around the extinction time problem, can't be the best way
NewSecondPeak = []
NewTimeSecondPeak = []

for j in range(trials):
    if TimeSecondPeak[j] > 200:
        NewSecondPeak.append(SecondPeak[j])
        NewTimeSecondPeak.append(TimeSecondPeak[j])

AvgSecondPeak = np.average(NewSecondPeak)
AvgTimeSecondPeak = np.average(NewTimeSecondPeak)

print("The second peak occurs (on average) at ", AvgTimeSecondPeak, " days, with the peak value: ", AvgSecondPeak)

# In[3]:

Quartile1 = np.percentile(SecondPeak, 25)
Quartile3 = np.percentile(SecondPeak, 75)
TQuartile1 = np.percentile(NewTimeSecondPeak, 25)
TQuartile3 = np.percentile(NewTimeSecondPeak, 75)

print("The IQR of the second fox peak TIME is between the range ", TQuartile1, " and ", TQuartile3, " days.")
print("The IQR of the second fox peak is between the range ", Quartile1, " and ", Quartile3, " foxes.")

# In[4]:

print("The probability that foxes die out is ", FDeathCount / trials * 100, "%.")