# In[1]:

# Kinetic Monte Carlo (KMC) Method

import numpy as np
from matplotlib import pyplot as plt
#import random

# Intializing Constants
k1 = 0.015       #days^-1 growth of rabbits
k2 = 0.00004     #days^-1foxes^-1 death of rabbits (eaten by foxes)
k3 = 0.0004      #days^-1rabbits^-1 growth of foxes (with eating rabbits)
k4 = 0.04        #days^-1 death of foxes
R0 = 400         #initial number of rabbits
F0 = 200         #initial number of foxes
days = 600       #to assign the end day

# Initializing Counter for Extinction Case
FDeathCount = 0     #will +1 when foxes die out in a trial
RDeathCount = 0     #will +1 when rabbits die out in a trial
trials = 1000       #total number of trials

# Initializing Lists for Second Peak (not "real") Calculations
SecondPeak = []         #to collect the 2nd peaks
TimeSecondPeak = []     #to collect the locations of the 2nd peak
             
for i in range(trials):
 
    # resetting below lists at the beginning of each trial
    KMCR = [R0]     # to collect the number of rabbits after a random event
    KMCF = [F0]     # to collect the number of foxes after the random event
    t = [0]         # to collect the time at which the random event happens
    
    while t[-1] < days:     #loop ends at the assigned end day

        # to calculate KMC parameters
        Rb = k1*KMCR[-1]                  #event of rabbit growth
        Rd = k2*KMCR[-1]*KMCF[-1]         #event of rabbit death
        Fb = k3*KMCR[-1]*KMCF[-1]         #event of fox growth
        Fd = k4*KMCF[-1]                  #event of fox death
        TotalRate = Rb + Rd + Fb + Fd     #sum of events

        Rando = np.random.rand()*TotalRate     #to generate the random number for deciding which event happens

        if KMCR[-1] == 0:     #to count in how many trials rabbits and foxes both die out
            RDeathCount = RDeathCount + 1
            # From Yujia: I add "FDeathCount = FDeathCount + 1" here, since as as long as rabbits die out, foxes must die out eventually.
            FDeathCount = FDeathCount + 1
            break
        elif KMCF[-1] == 0:     #to count in how many trials only foxes die out
            FDeathCount = FDeathCount + 1
            break

        # if neither rabbits nor foxes die out
        if Rando < Rb:     #the event of rabbit growth happens
            KMCR.append(KMCR[-1]+1)
            KMCF.append(KMCF[-1])
            t.append(t[-1]+((1/TotalRate)*np.log(1/np.random.rand())))
        elif Rando < (Rb + Rd):     #the event of rabbit death happens
            KMCR.append(KMCR[-1]-1)
            KMCF.append(KMCF[-1])
            t.append(t[-1]+((1/TotalRate)*np.log(1/np.random.rand())))
        elif Rando < (Rb + Rd + Fb):     #the event of fox growth happens
            KMCR.append(KMCR[-1])
            KMCF.append(KMCF[-1]+1)
            t.append(t[-1]+((1/TotalRate)*np.log(1/np.random.rand())))
        elif Rando < TotalRate:     #the event of fox death happens
            KMCR.append(KMCR[-1])
            KMCF.append(KMCF[-1]-1)
            t.append(t[-1]+((1/TotalRate)*np.log(1/np.random.rand())))

    # Getting the Second Peaks (not "Real") and Corresponding Locations for each Trial
    KMCSecondHalf = KMCF[len(KMCF)//2:]              #the 2nd half of fox population (has equal number of or one more point than the 1st half) for each trial
    TSecondHalf = t[len(KMCF)//2:]                   #the 2nd half of time (has equal number of or one more point than the 1st half) for each trial
    Peak2 = np.max(KMCSecondHalf)                    #to look for the 2nd peak for each trial
    SecondPeak.append(Peak2)                         #to append the value of the 2nd peak to the end of the collecting list for each trial
    Timeslot = KMCSecondHalf.index(Peak2)            #to look for the index of 2nd peak for each trial
    TimeSecondPeak.append(TSecondHalf[Timeslot])     #to look for and append the value of the 2nd peak location to the end of the collecting list for each trial

    # Plotting rabbit (in blue) and fox (in red) population as a function of time
    plt.plot(t,KMCR,'b')
    plt.plot(t,KMCF,'r')
#    plt.legend(loc='best')

# Adding Plot Elements
plt.xlabel('Days')
plt.ylabel('Population')
plt.legend(['Rabbit', 'Fox'])
plt.show()
print("Foxes died out ",FDeathCount," times and Rabbits died out ",RDeathCount," times." )


# In[2]:

# From Yujia: I think you are not collecting the "real" 2nd peak above. Based on the code, no matter whether either rabbits or foxes die out, it will
#             continue to split "KMCF" and "t" lists in half, and get the peak as well as the corresponding loaction. For example, if it runs 100 trials,
#             there will be 100 numbers in both "SecondPeak" and "TimeSecondPeak" lists. So for calculating the average of the 2nd peak, you should use
#             the similar approach as calculating the average time.

#AvgSecondPeak = np.average(SecondPeak)

# Getting the "Real" Second Peaks and Corresponding Locations
NewSecondPeak = []         #to create a list for collecting the "real" 2nd peaks
NewTimeSecondPeak = []     #to create a list for collecting the "real" 2nd peak locations

for j in range (trials):     #if the corresponding time < 200 days, the 2nd peak is not the "real"
    if TimeSecondPeak[j]>200:
        NewSecondPeak.append(SecondPeak[j])             #to append the value of the "real" 2nd peak to the end of the collecting list
        NewTimeSecondPeak.append(TimeSecondPeak[j])     #to append the value of the "real" 2nd peak location to the end of the collecting list

AvgSecondPeak = np.average(NewSecondPeak)             #to calculate the average of the "real" 2nd peaks
AvgTimeSecondPeak = np.average(NewTimeSecondPeak)     #to calculate the average of the "real" 2nd peak locations

# Outputs
print("The second peak occurs (on average) at ", AvgTimeSecondPeak, " days, with the peak value: ", AvgSecondPeak, ".")


# In[3]:

# Calculating Interquartile Ranges
Quartile1 = np.percentile(SecondPeak,25)
Quartile3 = np.percentile(SecondPeak,75)
TQuartile1 = np.percentile(NewTimeSecondPeak,25)
TQuartile3 = np.percentile(NewTimeSecondPeak,75)

# Outputs
print("The IQR of the second fox peak TIME is between the range ", TQuartile1, " and ", TQuartile3, " days.")
print("The IQR of the second fox peak is between the range ", Quartile1, " and ", Quartile3, " foxes.")


# In[4]:

# Outputs
print("The probability that foxes die out is ", float(FDeathCount)/float(trials)*100, "%.")     #has to convert "FDeathCount" and "trials" to float