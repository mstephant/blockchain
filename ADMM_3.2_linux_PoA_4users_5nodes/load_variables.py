#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 11:24:43 2020

@author: matthieu

- Load variables from result files 
- Plot

With pickle, variables should be loaded exactly in the same order that they have been saved. 
So all variables have to be loaded. 

"""
#%% Import variables

import pickle
import matplotlib.pyplot as plt
import pandas as pd


filename = 'Results_Agent'
with open(filename, 'rb') as f:
  timeVector = pickle.load(f) 
  price_imp = pickle.load(f)  
  coeff_prix_PV = pickle.load(f)
  coeff_prix_B = pickle.load(f)
  X_historic = pickle.load(f)
  L_historic = pickle.load(f)
  R_historic = pickle.load(f)
  S_historic = pickle.load(f)
  rho_historic = pickle.load(f)
  iterations = pickle.load(f)
  Pmax = pickle.load(f)
  gasSpend_agent = pickle.load(f)
  gasSpend_without_initialization_agent = pickle.load(f)
  etherSpend_agent = pickle.load(f)
X_historic_agent = X_historic
L_historic_agent = L_historic

length = len(timeVector)



    
    
#
#filename = 'Results_EV1'
#with open(filename, 'rb') as f:
#  X_historic_EV1 = pickle.load(f)
#  SOC_EV1 = pickle.load(f)
#  gasSpend_EV1 = pickle.load(f)
#  gasSpend_without_initialization_EV1 = pickle.load(f)
#  etherSpend_EV1 = pickle.load(f)
#  L_historic = pickle.load(f)
#L_historic_EV1 = L_historic
#X_EV1_final = pd.DataFrame({'p':-X_historic_EV1[-1][0:length],\
#                            'q':X_historic_EV1[-1][length:2*length],\
#                            '$\lambda$':X_historic_EV1[-1][2*length:3*length],\
#                            'r':X_historic_EV1[-1][3*length:4*length]},index = timeVector)
#SOC_EV1 = pd.DataFrame(SOC_EV1, index = timeVector)



filename = 'Results_EV2'
with open(filename, 'rb') as f:
  X_historic_EV2 = pickle.load(f)
  SOC_EV2 = pickle.load(f)
  gasSpend_EV2 = pickle.load(f)
  gasSpend_without_initialization_EV2 = pickle.load(f)
  etherSpend_EV2 = pickle.load(f)
  L_historic = pickle.load(f)
L_historic_EV2 = L_historic
X_EV2_final = pd.DataFrame({'p':-X_historic_EV2[-1][0:length],\
                            'q':X_historic_EV2[-1][length:2*length],\
                            '$\lambda$':X_historic_EV2[-1][2*length:3*length],\
                            'r':X_historic_EV2[-1][3*length:4*length]},index = timeVector)
SOC_EV2 = pd.DataFrame(SOC_EV2, index = timeVector)



#filename = 'Results_consumer_5RNS'
#with open(filename, 'rb') as f:
#  X_historic_consumer_5RNS = pickle.load(f)
#  forecast_consumer_5RNS = pickle.load(f)
#  gasSpend_consumer_5RNS = pickle.load(f) 
#  gasSpend_without_initialization_consumer_5RNS = pickle.load(f)
#  etherSpend_consumer_5RNS = pickle.load(f)
#  L_historic = pickle.load(f)
#L_historic_consumer_5RNS = L_historic
#X_consumer_5RNS_final = pd.DataFrame({'p':-X_historic_consumer_5RNS[-1][0:length], \
#                                      'q':X_historic_consumer_5RNS[-1][length:2*length],\
#                                      '$\lambda$':X_historic_consumer_5RNS[-1][2*length:3*length],\
#                                      'r':X_historic_consumer_5RNS[-1][3*length:4*length]},index = timeVector)


#filename = 'Results_consumer_HA'
#with open(filename, 'rb') as f:
#  X_historic_consumer_HA= pickle.load(f)
#  forecast_consumer_HA = pickle.load(f)
#  gasSpend_consumer_HA= pickle.load(f) 
#  etherSpend_consumer_HA = pickle.load(f)
#  gasSpend_without_initialization_consumer_HA = pickle.load(f)
#  L_historic = pickle.load(f)
#L_historic_consumer_HA = L_historic
#X_consumer_HA_final = pd.DataFrame({'p':-X_historic_consumer_HA[-1][0:length], \
#                                     'q':X_historic_consumer_HA[-1][length:2*length],\
#                                     '$\lambda$':X_historic_consumer_HA[-1][2*length:3*length],\
#                                     'r':X_historic_consumer_HA[-1][3*length:4*length],\
#                                     'forecast': forecast_consumer_HA},index = timeVector)
 

    
filename = 'Results_consumer_HEI'
with open(filename, 'rb') as f:
  X_historic_consumer_HEI= pickle.load(f)
  forecast_consumer_HEI = pickle.load(f)
  gasSpend_consumer_HEI= pickle.load(f) 
  gasSpend_without_initialization_consumer_HEI = pickle.load(f)
  etherSpend_consumer_HEI = pickle.load(f)
  L_historic = pickle.load(f)
L_historic_consumer_HEI = L_historic
X_consumer_HEI_final = pd.DataFrame({'p':-X_historic_consumer_HEI[-1][0:length], \
                                     'q':X_historic_consumer_HEI[-1][length:2*length],\
                                     '$\lambda$':X_historic_consumer_HEI[-1][2*length:3*length],\
                                     'r':X_historic_consumer_HEI[-1][3*length:4*length],\
                                     'forecast': forecast_consumer_HEI},index = timeVector)
 
    
#filename = 'Results_consumer_RIZOMME'
#with open(filename, 'rb') as f:
#  X_historic_consumer_RIZOMME = pickle.load(f)
#  forecast_consumer_RIZOMME = pickle.load(f)
#  gasSpend_consumer_RIZOMME= pickle.load(f) 
#   gasSpend_without_initialization_consumer_RIZOMME = pickle.load(f)
#  etherSpend_consumer_RIZOMME = pickle.load(f)
#  L_historic = pickle.load(f)
#L_historic_consumer_RIZOMME = L_historic
#X_consumer_RIZOMME_final = pd.DataFrame({'p':-X_historic_consumer_RIZOMME[-1][0:length], \
#                                     'q':X_historic_consumer_RIZOMME[-1][length:2*length],\
#                                     '$\lambda$':X_historic_consumer_RIZOMME[-1][2*length:3*length],\
#                                     'r':X_historic_consumer_RIZOMME[-1][3*length:4*length],\
#                                     'forecast': forecast_consumer_RIZOMME},index = timeVector)

#
#filename = 'Results_PV_5RNS'
#with open(filename, 'rb') as f:
#    X_historic_PV_5RNS = pickle.load(f)
#    forecast_PV_5RNS = pickle.load(f)
#    gasSpend_PV_5RNS = pickle.load(f)
#    gasSpend_without_initialization_PV_5RNS = pickle.load(f)
#    etherSpend_PV_5RNS = pickle.load(f)
#    L_historic = pickle.load(f)
#L_historic_PV_5RNS = L_historic
#X_consumer_PV_5RNS_final = pd.DataFrame({'p':X_historic_consumer_PV_5RNS[-1][0:length], \
#                                     '-q':-X_historic_consumer_PV_5RNS[-1][length:2*length],\
#                                     '$\lambda$':X_historic_consumer_PV_5RNS[-1][2*length:3*length],\
#                                     '-r':-X_historic_consumer_PV_5RNS[-1][3*length:4*length],\
#                                     'forecast': forecast_consumer_PV_5RNS},index = timeVector)


filename = 'Results_PV_RIZOMME'
with open(filename, 'rb') as f:
    X_historic_PV_RIZOMME = pickle.load(f)
    forecast_PV_RIZOMME = pickle.load(f)
    gasSpend_PV_RIZOMME = pickle.load(f)
    gasSpend_without_initialization_PV_RIZOMME = pickle.load(f)
    etherSpend_PV_RIZOMME = pickle.load(f)
    L_historic = pickle.load(f)
L_historic_PV_RIZOMME = L_historic
X_PV_RIZOMME_final = pd.DataFrame({'p':X_historic_PV_RIZOMME[-1][0:length], \
                                     '-q':-X_historic_PV_RIZOMME[-1][length:2*length],\
                                     '$\lambda$':X_historic_PV_RIZOMME[-1][2*length:3*length],\
                                     '-r': -X_historic_PV_RIZOMME[-1][3*length:4*length],\
                                     'forecast':   forecast_PV_RIZOMME},index = timeVector)




filename = 'Results_storage'
with open(filename, 'rb') as f:
  X_historic_storage = pickle.load(f)
  SOC_storage = pickle.load(f)
  gasSpend_Storage = pickle.load(f)
  gasSpend_without_initialization_Storage = pickle.load(f)
  etherSpend_Storage = pickle.load(f)
  L_historic = pickle.load(f)
L_historic_storage = L_historic
X_storage_final = pd.DataFrame({'-p':-X_historic_storage[-1][0:length],\
                                'q':X_historic_storage[-1][length:2*length],\
                                '$\lambda$':X_historic_storage[-1][2*length:3*length],\
                                'r':X_historic_storage[-1][3*length:4*length]},index = timeVector)
SOC_storage = pd.DataFrame(SOC_storage, index = timeVector)


#%% Plot



# EV2
plt.figure()
X_EV2_final.plot(ylabel = "Power (kW)", title="Profile of EV$_2$")
plt.savefig('user_EV2_power.pdf')  

plt.figure()
SOC_EV2.plot(ylabel = "SOC (%)", title="SOC EV$_2$")
plt.savefig('user_EV2_SOC.pdf')  

# consumer_HEI
plt.figure()
X_consumer_HEI_final.plot(ylabel = "Power (kW)", title="Profile of building HEI", style = ['-','-','-','-','--'])
plt.savefig('user_consumer_HEI.pdf')  

# PV_RIZOMME
plt.figure()
X_PV_RIZOMME_final.plot(ylabel = "Power (kW)", title="Profile of PV 5RNS", style = ['-','-','-','-','--'])
plt.savefig('user_PV_RIZOMME.pdf')  

# storage
plt.figure()
X_storage_final.plot(ylabel = "Power (kW)", title="Profile of battery")
plt.savefig('user_storage_power.pdf')  

plt.figure()
SOC_storage.plot(ylabel = "SOC (%)", title="SOC battery")
plt.savefig('user_storage_SOC.pdf')  

#%% gas spend
gasSpend_EV2_per_iteration = gasSpend_EV2/iterations
gasSpend_PV_RIZOMME_per_iteration = gasSpend_PV_RIZOMME/iterations
gasSpend_Storage_per_iteration = gasSpend_Storage/iterations
gasSpend_consumer_HEI_per_iteration = gasSpend_consumer_HEI/iterations
gasSpend_agent_per_iteration = gasSpend_agent/iterations


gasSpend_moyen_per_iteration = (gasSpend_EV2_per_iteration + gasSpend_PV_RIZOMME_per_iteration + gasSpend_Storage_per_iteration  + gasSpend_consumer_HEI_per_iteration) /4


# without initialization
gasSpend_without_initialization_EV2_per_iteration = gasSpend_without_initialization_EV2/iterations
gasSpend_without_initialization_PV_RIZOMME_per_iteration = gasSpend_without_initialization_PV_RIZOMME/iterations
gasSpend_without_initialization_Storage_per_iteration = gasSpend_without_initialization_Storage/iterations
gasSpend_without_initialization_consumer_HEI_per_iteration = gasSpend_without_initialization_consumer_HEI/iterations
gasSpend_without_initialization_agent_per_iteration = gasSpend_without_initialization_agent/iterations


gasSpend_moyen_per_iteration = (gasSpend_EV2_per_iteration + gasSpend_PV_RIZOMME_per_iteration + gasSpend_Storage_per_iteration  + gasSpend_consumer_HEI_per_iteration) /4






