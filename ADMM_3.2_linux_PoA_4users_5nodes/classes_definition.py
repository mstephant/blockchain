# -*- coding: utf-8 -*-
"""
Last modification: 15/09/2020
@author: matthieu.stephant

Definition of users classes (EV, storage, PV generator, load)
"""

import utility
import numpy as np
import pandas as pd


###############################################################################

class User:
    """General user class:
            - userType (1:EV, 2: storage, 3: PV generator, 4: load)
            - optimizedProfile (np.array) """
            
    def __init__(self, userType):
        """Constructor """
        self.userType = userType
        self.optimizedProfile = np.array([])

    def updateProfile(self, newProfile):
        self.optimizedProfile = newProfile
        
###############################################################################    
    
class EV(User):
    """Class for EV, inheriting from User:
        - alpha (preferences) = [costs, local exchange, comfort]
        - chargingTime = [arrivalTime, departureTime]
        - Parameters: [PmaxEV,EmaxEV,etaEV,SOCmin,SOCmax,SOCinit]
        - SOCrequired 
        - SOC (np.array)"""
    
    def __init__(self):
        User.__init__(self,1) # constructor of User, with userType = 1 for EV  
            
    def setPreferences(self, costs, exchanges, comfort):
        self.alpha = [costs, exchanges, comfort]
        
        for i in self.alpha:
            if i<0 or i>1:
                raise ValueError('Error, each coefficient should be positive and <1.')
        
        if sum(self.alpha)!=1:
            raise ValueError('Error, sum of coefficients should be egal to 1.')
            
    def setParameters(self, PmaxEV,EmaxEV,etaEV,SOCmin,SOCmax,SOCinit):
        """
        PmaxEV (kW), EmaxEV (kWh)
        """
        self.parameters = [PmaxEV,EmaxEV,etaEV,SOCmin,SOCmax,SOCinit]
        
    def setChargingTime(self, arrivalTime, departureTime):
        if arrivalTime > departureTime:
            raise ValueError('Error, arrivalTime should be less than departureTime')
        self.arrivalTime = arrivalTime
        self.departureTime = departureTime
        self.chargingTime = departureTime - arrivalTime
        self.chargingTime = self.chargingTime.total_seconds()/3600 # chargingTime in hours
        
    def setSOCrequired(self,SOCrequired):
        # maximum SOC possible at the end of charging time
        SOClimit =  self.parameters[5]+(self.parameters[2]*self.parameters[0]* \
        self.chargingTime)/self.parameters[1]
        # SOClimit = SOCinit + eta*Pmax*chargingTime/Emax
        
        if SOCrequired>self.parameters[4]:
            raise ValueError('Error: SOCrequired > SOCmax')

        elif SOCrequired>SOClimit:
            raise ValueError('Error: SOCrequired_EV1 too high. SOC cannot reach this level during the time of charge')       
        self.SOCrequired = SOCrequired
        
    def optimization(self,startTime, endTime, timestep, price_imp, price_PV, price_B, rho, L, yPVforecast):
        """ Optimization for a specific time interval.
            Use the 'optimEV' function 
            Inputs:
                - startTime, endTime: pd.to_datetime()
                - timestep: minutes
                - price_imp, price_PV, price_B, yPVforecast: pd.Series
                - rho: floats
                - L : np.array
        """
 
        newProfile = utility.optimEV(self.alpha,self.arrivalTime,
                                          self.departureTime,self.parameters,
                                          self.SOCrequired, startTime, endTime, 
                                          timestep, price_imp, price_PV, price_B, rho, L, 
                                          self.optimizedProfile, yPVforecast)
        self.updateProfile(newProfile)
        
        
    def CalculateSOC(self,timestep):
        """ Calculation of SOC, with the profile saved into the class.
               - timestep should be in hours !!!!!
        """
        length = int(len(self.optimizedProfile)/4)
        SOC = np.zeros(length)
        SOC[0] = self.parameters[5]
        for i in range(1,length):
            SOC[i]=SOC[i-1]-self.parameters[2]*(self.optimizedProfile[i-1])*timestep/self.parameters[1] 
        self.SOC = SOC
            

        
###############################################################################     
        
        
class Storage(User):
    """Class for storage, inheriting from User:
        - beta (preferences) = [cost, local exchanges]
        - Parameters: [PmaxIn_b,PmaxOut_b,Emax_b,eta_b,SOCmin_b,SOCmax_b,SOCinit_b]
        - SOC (np.array)
    
    """
    
    def __init__(self):
        User.__init__(self,2) # constructor of User, with userType = 2 for Storage 
            
    def setPreferences(self, costs, exchanges):
        self.beta = [costs, exchanges]
        
        for i in self.beta:
            if i<0 or i>1:
                raise ValueError('Error, each coefficient should be positive and <1.')
        
        if sum(self.beta)!=1:
            raise ValueError('Error, sum of coefficients should be egal to 1.')
            
    def setParameters(self, PmaxIn_b,PmaxOut_b,Emax_b,eta_b,SOCmin_b,SOCmax_b,SOCinit_b):
        self.parameters = [PmaxIn_b,PmaxOut_b,Emax_b,eta_b,SOCmin_b,SOCmax_b,SOCinit_b]
        
    def optimization(self,startTime, endTime, timestep,  price_imp, price_PV, price_B, rho, L, yPVforecast):
        """ Optimization for a specific time interval.
            Use the 'optimSt' function 
            Inputs:
                - startTime, endTime: pd.to_datetime()
                - timestep: minutes
                - price_imp, price_PV,  yPVforecast : pd.Series
                - L, rho: floats
        """
        newProfile = utility.optimStor(self.beta,self.parameters,startTime,
                                     endTime,timestep,price_imp,price_PV,price_B,rho,L,
                                     self.optimizedProfile, yPVforecast)
        
        self.updateProfile(newProfile)
        
    def CalculateSOC(self,timestep):
        """ Calculation of SOC, with the profile saved into the class.
        - timestep should be in hours
        """
        length = int(len(self.optimizedProfile)/4)
        SOC = np.zeros(length)
        SOC[0] = self.parameters[6]
        for i in range(1,length):
            SOC[i]=SOC[i-1]-self.parameters[3]*self.optimizedProfile[i-1]*timestep/self.parameters[2]
        self.SOC = SOC
     
         
###############################################################################        

class GeneratorPV(User):
    """Class for PV generator, inheriting from User:
        - gamma (preferences) = [cost, comfort]"""
    
    def __init__(self):
        User.__init__(self,3) # constructor of User, with userType = 3 for GeneratorPV  
            
    def setPreferences(self, costs, exchanges,comfort):
        self.gamma = [costs, exchanges,comfort]
        
        for i in self.gamma:
            if i<0 or i>1:
                raise ValueError('Error, each coefficient should be positive and <1.')
        
        if sum(self.gamma)!=1:
            raise ValueError('Error, sum of coefficients should be egal to 1.')

    def setForecast(self, forecastedProfile):
        """
        forecastedProfile: pd.Series
        """
        self.forecast = forecastedProfile
    
    def optimization(self,startTime, endTime, timestep, price_imp, price_PV, rho, L):
        """ Optimization for a specific time interval.
            Use the 'optimPV' function 
            Inputs:
                - startTime, endTime: pd.to_datetime()
                - timestep: minutes
                - price_imp, price_PV: pd.Series
                - L, rho: floats
        """
        
        newProfile = utility.optimPV(self.gamma,startTime,endTime,timestep,
                                     self.forecast,price_imp, price_PV, rho,L,
                                     self.optimizedProfile)
        self.updateProfile(newProfile)
                
###############################################################################

class Consumer(User): 
    """Class for consumer, inheriting from User:
        - delta (preferences) = (costs, local exchanges,comfort)
        - flexibility
    """
    
    def __init__(self):
        User.__init__(self,4) # constructor of User, with userType = 4 for Load
            
    def setPreferences(self, costs, exchanges, comfort):
        self.delta = [costs, exchanges, comfort]
        
        for i in self.delta:
            if i<0 or i>1:
                raise ValueError('Error, each coefficient should be positive and <1.')
        
        if sum(self.delta)!=1:
            raise ValueError('Error, sum of coefficients should be egal to 1.')
            
    def setForecast(self, forecastedProfile):
        """
        forecastedProfile: pd.Series
        """
        self.forecast = forecastedProfile
        
    def setFlexibility(self, flexibility):
        self.flexibility = flexibility    
   

    def optimization(self,startTime, endTime, timestep, price_PV, price_imp, price_B, rho, L, yPVforecast):
        """ Optimization for a specific time interval.
            Use the 'optimConsumer' function 
            Inputs:
                - startTime, endTime: pd.to_datetime()
                - timestep: minutes
                - price_PV, price_imp,  yPVforecast : pd.Series
                - L, rho: floats
        """
        newProfile = utility.optimConsumer(self.delta,self.flexibility,self.forecast,
                                       startTime,endTime,timestep,price_PV,price_imp,price_B,
                                       rho,L,self.optimizedProfile, yPVforecast)

        self.updateProfile(newProfile)
                