
# coding: utf-8

# # Rabbits and foxes
# 
# There are initially 400 rabbits and 200 foxes on a farm (but it could be two cell types in a 96 well plate or something, if you prefer bio-engineering analogies). Plot the concentration of foxes and rabbits as a function of time for a period of up to 600 days. The predator-prey relationships are given by the following set of coupled ordinary differential equations:
# 
# \begin{align}
# \frac{dR}{dt} &= k_1 R - k_2 R F \tag{1}\\
# \frac{dF}{dt} &= k_3 R F - k_4 F \tag{2}\\
# \end{align}
# 
# * Constant for growth of rabbits $k_1 = 0.015$ day<sup>-1</sup>
# * Constant for death of rabbits being eaten by foxes $k_2 = 0.00004$ day<sup>-1</sup> foxes<sup>-1</sup>
# * Constant for growth of foxes after eating rabbits $k_3 = 0.0004$ day<sup>-1</sup> rabbits<sup>-1</sup>
# * Constant for death of foxes $k_1 = 0.04$ day<sup>-1</sup>
# 
# Also plot the number of foxes versus the number of rabbits.
# 
# Then try also with 
# * $k_3 = 0.00004$ day<sup>-1</sup> rabbits<sup>-1</sup>
# * $t_{final} = 800$ days
# 
# *This problem is based on one from Chapter 1 of H. Scott Fogler's textbook "Essentials of Chemical Reaction Engineering".*
# 

# # Solving ODEs
# 
# *Much of the following content reused under Creative Commons Attribution license CC-BY 4.0, code under MIT license (c)2014 L.A. Barba, G.F. Forsyth. Partly based on David Ketcheson's pendulum lesson, also under CC-BY. https://github.com/numerical-mooc/numerical-mooc*
# 
# Let's step back for a moment. Suppose we have a first-order ODE $u'=f(u)$. You know that if we were to integrate this, there would be an arbitrary constant of integration. To find its value, we do need to know one point on the curve $(t, u)$. When the derivative in the ODE is with respect to time, we call that point the _initial value_ and write something like this:
# 
# $$u(t=0)=u_0$$
# 
# In the case of a second-order ODE, we already saw how to write it as a system of first-order ODEs, and we would need an initial value for each equation: two conditions are needed to determine our constants of integration. The same applies for higher-order ODEs: if it is of order $n$, we can write it as $n$ first-order equations, and we need $n$ known values. If we have that data, we call the problem an _initial value problem_.
# 
# Remember the definition of a derivative? The derivative represents the slope of the tangent at a point of the curve $u=u(t)$, and the definition of the derivative $u'$ for a function is:
# 
# $$u'(t) = \lim_{\Delta t\rightarrow 0} \frac{u(t+\Delta t)-u(t)}{\Delta t}$$
# 
# If the step $\Delta t$ is already very small, we can _approximate_ the derivative by dropping the limit. We can write:
# 
# $$\begin{equation}
# u(t+\Delta t) \approx u(t) + u'(t) \Delta t
# \end{equation}$$
# 
# With this equation, and because we know $u'(t)=f(u)$, if we have an initial value, we can step by $\Delta t$ and find the value of $u(t+\Delta t)$, then we can take this value, and find $u(t+2\Delta t)$, and so on: we say that we _step in time_, numerically finding the solution $u(t)$ for a range of values: $t_1, t_2, t_3 \cdots$, each separated by $\Delta t$. The numerical solution of the ODE is simply the table of values $t_i, u_i$ that results from this process.
# 

# # Euler's method
# *Also known as "Simple Euler" or sometimes "Simple Error".*
# 
# The approximate solution at time $t_n$ is $u_n$, and the numerical solution of the differential equation consists of computing a sequence of approximate solutions by the following formula, based on Equation (10):
# 
# $$u_{n+1} = u_n + \Delta t \,f(u_n).$$
# 
# This formula is called **Euler's method**.
# 
# For the equations of the rabbits and foxes, Euler's method gives the following algorithm that we need to implement in code:
# 
# \begin{align}
# R_{n+1} & = R_n + \Delta t \left(k_1 R_n - k_2 R_n F_n \right) \\
# F_{n+1} & = F_n + \Delta t \left( k_3 R_n F_n - k_4 F_n \right).
# \end{align}
# 

# In[1]:

get_ipython().magic('matplotlib inline')
import numpy as np
from matplotlib import pyplot as plt


#Intializing Constants
k1 = 0.015 #days^-1 Growth of Rabbits
k2 = 0.00004 #days^-1foxes^-1 Rabbits eaten by foxes
k3 = 0.0004 #days^-1 foxes eating rabbits
k4 = 0.04  #days^-1 Death of foxes
R0 = 400 #Initial Rabbits
F0 = 200 #Initial Foxes
days = 600
time = np.linspace(0,days,600)


# In[2]:

#Euler's Method

Rabbits = [R0]
Foxes = [F0]
dt = 1 
i = 0

for i in range(0,599):

    Rabbits.append(Rabbits[i] + dt*(k1*Rabbits[i] - k2*Rabbits[i]*Foxes[i]))
    Foxes.append(Foxes[i] + dt*(k3*Rabbits[i]*Foxes[i] - k4*Foxes[i]))
  
plt.plot(time,Rabbits)
plt.plot(time,Foxes)
plt.show()


# In[3]:

#ODEint Method

from scipy.integrate import odeint

def ODE(RF,t,k1,k2,k3,k4):
    R,F=RF
    dRdt = k1*R - k2*R*F
    dFdt = k3*R*F-k4*F
    bothDiff = [dRdt,dFdt]
    return bothDiff

RF0 = [R0,F0]
ans = odeint(ODE,RF0,time,args=(k1,k2,k3,k4))

plt.plot(time,ans[:,0],label='Rabbits') #Not Sure how this array splitting works
plt.plot(time,ans[:,1],label='Foxes') #Not sure how this array splitting works
plt.legend(loc='best')
plt.show()


# In[4]:

#Finding Second Max from ODEint

ODERabbits = ans[:,0]
ODEFoxes = ans[:,1]

FirstHalf = ODEFoxes[:len(ODEFoxes)//2]
SecondHalf = ODEFoxes[len(ODEFoxes)//2:]
MaxFoxes1 = np.max(FirstHalf)
MaxFoxes2 = np.max(SecondHalf)
TimeMaxFoxes1 = np.argmax(FirstHalf)
TimeMaxFoxes2 = len(ODEFoxes)//2+np.argmax(SecondHalf)


print("The first max: ",MaxFoxes1," occurs at ",TimeMaxFoxes1," days.")
print("The second max: ",MaxFoxes2," occurs at ",TimeMaxFoxes2," days.")


# In[22]:

#Kinetic Monte Carlo Method

import random

#Initialize Counter for Extinction Case
FDeathCount = 0
RDeathCount = 0
trials = 1000

#Initialize Arrays for Second Peak Calculations
SecondPeak = []
TimeSecondPeak = []
             
for i in range(trials):
 
    #Reset for Each Trial

    t= [0]
    KMCR = [R0]
    KMCF = [F0]
    
    while t[-1] < days:


        #Rates multiply by the last item in the array (Most Recent Population)
        Rb = k1*KMCR[-1]
        Rd = k2*KMCR[-1]*KMCF[-1]
        Fb = k3*KMCR[-1]*KMCF[-1]
        Fd = k4*KMCF[-1]
        TotalRate = Rb + Rd + Fb + Fd

        #Generate Random Number
        Rando = np.random.rand()*TotalRate


        if KMCR[-1] == 0:
            #Rabbits are Dead
            RDeathCount = RDeathCount + 1
            break
        elif KMCF[-1] == 0:
            #Foxes are Dead
            FDeathCount = FDeathCount + 1
            break
        if Rando < Rb:
            #Rabbit is Born
            KMCR.append(KMCR[-1]+1)
            KMCF.append(KMCF[-1])
            t.append(t[-1]+((1/TotalRate)*np.log(1/np.random.rand())))
        elif Rando < (Rb + Rd):
            #Rabbit Dies
            KMCR.append(KMCR[-1]-1)
            KMCF.append(KMCF[-1])
            t.append(t[-1]+((1/TotalRate)*np.log(1/np.random.rand())))
        elif Rando < (Rb + Rd + Fb):
            #Fox is Born
            KMCR.append(KMCR[-1])
            KMCF.append(KMCF[-1]+1)
            t.append(t[-1]+((1/TotalRate)*np.log(1/np.random.rand())))
        elif Rando < TotalRate:
            #FoxDies
            KMCR.append(KMCR[-1])
            KMCF.append(KMCF[-1]-1)
            t.append(t[-1]+((1/TotalRate)*np.log(1/np.random.rand())))

    
    #Getting Second Peaks and Second Peak Times for each Trial
    KMCSecondHalf = KMCF[len(KMCF)//2:]
    TSecondHalf = t[len(KMCF)//2:]     
    Peak2 = np.max(KMCSecondHalf)
    SecondPeak.append(Peak2)
    Timeslot = KMCSecondHalf.index(Peak2)
    TimeSecondPeak.append(TSecondHalf[Timeslot])

    #Stole plotting on the same graph from @Justin
    plt.plot(t,KMCR,'b')
    plt.plot(t,KMCF,'r')
    plt.legend(loc='best')

plt.xlabel('Days')
plt.ylabel('Population')
plt.show()
print("Foxes died out ",FDeathCount," times and Rabbits died out ",RDeathCount," times." )


# In[23]:

AvgSecondPeak = np.average(SecondPeak)

#My way of getting around the extinction time problem, can't be the best way
NewTimeSecondPeak = []
for j in range (trials-1):
    if TimeSecondPeak[j]>200:
        NewTimeSecondPeak.append(TimeSecondPeak[j])
AvgTimeSecondPeak = np.average(NewTimeSecondPeak)

print("The second peak occurs (on average) at ",AvgTimeSecondPeak," days, with the peak value: ",AvgSecondPeak)


# In[24]:

Quartile1 = np.percentile(SecondPeak,25)
Quartile3 = np.percentile(SecondPeak,75)
TQuartile1 = np.percentile(NewTimeSecondPeak,25)
TQuartile3 = np.percentile(NewTimeSecondPeak,75)

print("The IQR of the second fox peak TIME is between the range ",TQuartile1," and ",TQuartile3," days.")
print("The IQR of the second fox peak is between the range ",Quartile1," and ",Quartile3," foxes.")


# In[25]:

print("The probability that foxes die out is ",FDeathCount/trials*100,"%.")

# pull request practice
print("Pull Request Practice")


# In[ ]:

#Lessons Learned

#Trial and error is your best friend. I've learned to read error messages fairly well.
#Seeing as I mess up quite a bit, Python tries very hard to tell you what to do.

#I had to search how to plot things on the same graph (and ended up looking at the easy solution on Justin's code)

#Numpy can do SO MANY THINGS, one day I'll use some of them.

#Kinetic Monte Carlo sounds really hard on the Wiki page, but it's actually very simple in theory.


#I still don't really understand ODEint, more specifically the input/output arguments.

