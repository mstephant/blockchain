# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 11:14:55 2020

@author: matthieu.stephant

Calculations of indicators for final results of distributed ADMM smart contract
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import pickle

### Loading variables from results files

filename = 'Results\Agent'
with open(filename, 'rb') as f:
    timeVector = pickle.load(f)
    price_imp = pickle.load(f)
    coeff_prix = pickle.load(f)
    X_historic = pickle.load(f)
    L_historic = pickle.load(f)
    R_historic = pickle.load(f)
    S_historic = pickle.load(f)
    rho_historic = pickle.load(f)
    iterations = pickle.load(f)
    Pmax = pickle.load(f)


filename = 'Results/consumer_5RNS'
with open(filename, 'rb') as f:
    X_historic_consumer_5RNS = pickle.load(f)
    forecast_consumer_5RNS = pickle.load(f)
    X_consumer_5RNS = pickle.load(f)

filename = 'Results/consumer_HA'
with open(filename, 'rb') as f:
    X_historic_consumer_HA = pickle.load(f)
    forecast_consumer_HA = pickle.load(f)
    X_consumer_HA = pickle.load(f)
    
filename = 'Results/consumer_HEI'
with open(filename, 'rb') as f:
    X_historic_consumer_HEI = pickle.load(f)
    forecast_consumer_HEI = pickle.load(f)
    X_consumer_HEI = pickle.load(f)
    
filename = 'Results/consumer_RIZOMME'
with open(filename, 'rb') as f:
    X_historic_consumer_RIZOMME = pickle.load(f)
    forecast_consumer_RIZOMME = pickle.load(f)
    X_consumer_RIZOMME = pickle.load(f)
    
    
filename = 'Results/PV_5RNS'
with open(filename, 'rb') as f:
    X_historic_PV_5RNS = pickle.load(f)
    forecast_PV_5RNS = pickle.load(f)
    X_PV_5RNS = pickle.load(f)

filename = 'Results/PV_RIZOMME'
with open(filename, 'rb') as f:
    X_historic_PV_RIZOMME = pickle.load(f)
    forecast_PV_RIZOMME = pickle.load(f)
    X_PV_RIZOMME = pickle.load(f)


filename = 'Results\EV1'
with open(filename, 'rb') as f:
    X_historic_EV1 = pickle.load(f)
    SOC_EV1 = pickle.load(f)
    X_EV1 = pickle.load(f)
    
filename = 'Results\EV2'
with open(filename, 'rb') as f:
    X_historic_EV2 = pickle.load(f)
    SOC_EV2 = pickle.load(f)
    X_EV2 = pickle.load(f)    
    
filename = 'Results\storage'
with open(filename, 'rb') as f:
    X_historic_storage = pickle.load(f)
    SOC_storage = pickle.load(f)
    X_SOC = pickle.load(f)

