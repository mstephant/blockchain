# -*- coding: utf-8 -*-
"""
Last modification: 22/06/2020
@author: matthieu.stephant

"""

import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl


#matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.rcParams['ps.fonttype'] = 42

# load the workspace:
with open('Results\Agent', 'rb') as f:
  X_historic =  pickle.load(f)

# plot
Pmax = 300 
constraint = Pmax*np.ones(73)


startTime = pd.to_datetime('2019-06-13 00:00')
endTime = pd.to_datetime('2019-06-14 00:00')
timeStep = 20 # minutes
timeVector = pd.date_range(startTime, endTime,freq = str(timeStep)+'T')
length = len(timeVector)


#fig, ax = plt.subplots()
#
#mpl.rcParams['text.usetex'] = True
#
#ax.plot(timeVector, constraint, 'r--', label='Global constraint')
#ax.plot(timeVector,X_historic[0]/10**6, 'b-', label='First iteration (no penalty)')
#ax.plot(timeVector,X_historic[-1]/10**6, 'g-', label='Final state')    
#legend = ax.legend(loc='upper right', fontsize='small')
#plt.xlabel('Date')
#plt.ylabel('Power (kW)')  
#plt.show()
#plt.savefig('Results\global_figure.pdf')  


# plot
fig,ax = plt.subplots()
plt.title('Power profile of building HEI')
ax.plot(X_consumer_HEI,label='Consumption', linestyle='-', color='b')
ax.plot(forecast_consumer_HEI,label='forecast_consumer_HEI',linestyle='--', color='g')
ax.set_xlabel('Date')
ax.set_ylabel('Power (kW)')
ax.legend(loc='upper right', fontsize='small')
plt.savefig('Results\consumer_HEI.pdf')  