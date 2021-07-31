# -*- coding: utf-8 -*-

"""
Last modification: 14/10/2020
@author: matthieu.stephant

script for storage 
interaction with blockchain node
PoW 
smart contract: ADMM_version3.2.sol
	-> calculation of newL is done in the python code and not in the smart contract


use of MATALB function util_func_B_3.m for solving optimization problem, because
does not converge with Python
"""

import json
import classes_definition as cl
import numpy as np
import pandas as pd
import tools as tl 
import pickle
#from web3.middleware import geth_poa_middleware
import time
import matlab.engine

# start MATLAB engine
eng = matlab.engine.start_matlab()

#from web3.auto import web3 # IPC automatic connection
from web3 import Web3
node_url = "http://localhost:8003"
web3 = Web3(Web3.HTTPProvider(node_url,request_kwargs={'timeout':500}))

# set pre-funded account as sender
web3.eth.defaultAccount = web3.eth.accounts[0]
web3.geth.personal.unlockAccount(web3.eth.defaultAccount,"0000",0) # unlock account permanently

# mining of all sealer nodes necessary
#web3.geth.miner.start()

# PoA compatibility
#web3.middleware_onion.inject(geth_poa_middleware, layer=0)

time.sleep(30)

# import contract details
with open("../../smart_contracts/ADMM_version3.2.json",'r') as f:
    abi = json.load(f)    
with open("address.txt", "r") as f:
    contract_address = f.read() 
    
contractADMM = web3.eth.contract(
    address=contract_address,
    abi=abi,
)

# instanciate class of user (EV, Storage, GeneratorPV, Consumer)
Storage = cl.Storage()

# gas counter initialization
gasSpend = 0
etherSpend = web3.eth.getBalance(web3.eth.defaultAccount)

# create a User struct
tx_hash = contractADMM.functions.newUser(Storage.userType).transact()
tx_receipt = web3.eth.waitForTransactionReceipt(tx_hash)
gasSpend += tx_receipt.gasUsed

# preferences
Storage.setPreferences(0.5,0.5)  # (cost, local exchanges)
    
# parameters
PmaxIn_b = 40
PmaxOut_b = 80
Emax_b = 203
eta_b = 97/100
SOCmin_b = 10/100
SOCmax_b = 100/100
SOCinit_b = 35/100
Storage.setParameters(PmaxIn_b,PmaxOut_b,Emax_b,eta_b,SOCmin_b,SOCmax_b,SOCinit_b)


# listen to event to know when to start optimization
l = 0
while l!=1:
    event_filter = contractADMM.events.StartOptim.createFilter(fromBlock=0)
    eventlist = event_filter.get_all_entries()
    l=len(eventlist)
print('Starting optimization')
    
# get time interval
[startTime, endTime, timeStep] = contractADMM.functions.getTimeInterval().call()
    # startTime and endTime are Unix time and timeStep in seconds
startTime = pd.to_datetime(startTime,unit='s')
endTime = pd.to_datetime(endTime,unit='s')
timeStep = timeStep/60 # from seconds to minutes
timeVector = pd.date_range(startTime, endTime,freq = str(timeStep)+'T')
length = len(timeVector)

# get price vector 
price_imp = contractADMM.functions.getPriceImp().call() 
price_imp = np.array(price_imp)
price_imp = pd.Series(price_imp/10**6,timeVector)

[coeff_pricePV, coeff_priceB] = contractADMM.functions.getPriceLocal().call() 
coeff_pricePV = coeff_pricePV/10**6
price_PV = pd.Series(data=coeff_pricePV*np.ones(length),index=timeVector)
coeff_priceB = coeff_priceB/10**6
price_B = pd.Series(data=coeff_priceB*np.ones(length),index=timeVector)


# listen to event to know when total PV forecast is available
l = 0
while l!=1:
    event_filter = contractADMM.events.SendForecast.createFilter(fromBlock=0)
    eventlist = event_filter.get_all_entries()
    l=len(eventlist)
print('Get total PV forecast')

# get PV forecast
yPVforecast = contractADMM.functions.getPVforecast().call()
yPVforecast = np.array(yPVforecast)
yPVforecast = pd.Series(yPVforecast/10**6,timeVector)


# initialisation
result = 0 
iterations = 0
X_historic = []
L_historic=[]

# re-define parameters for good MATLAB type
beta1 = eng.double(Storage.beta[0])
beta2 = eng.double(Storage.beta[1])
param = np.array(Storage.parameters)
parametersB = eng.cell2mat(param.tolist())
timestepM = eng.minutes(timeStep, nargout=1)
startTimeM = eng.datetime(startTime.year,startTime.month,startTime.day,startTime.hour,startTime.minute,startTime.second,nargout=1)
endTimeM = eng.datetime(endTime.year,endTime.month,endTime.day,endTime.hour,endTime.minute,endTime.second,nargout=1)
timeVectorM = eng.colon(startTimeM,timestepM,endTimeM,nargout=1)
timeVectorM = eng.transpose(timeVectorM)
price_PV_M = eng.timetable(timeVectorM,eng.transpose(eng.cell2mat(price_PV.values.tolist())),nargout=1)
price_b_M = eng.timetable(timeVectorM,eng.transpose(eng.cell2mat(price_B.values.tolist())),nargout=1)
price_imp_M = eng.timetable(timeVectorM,eng.transpose(eng.cell2mat(price_imp.values.tolist())),nargout=1)
p_forecast_PV_total_M = eng.transpose(eng.cell2mat(yPVforecast.values.tolist()))

gasSpend_initialization = gasSpend

while result ==0: # while stopping criteria of global optimization not reached
    
    time.sleep(15)###
    
    # get L and rho
    [L,rho,curtailPV] = contractADMM.functions.getOptimizationParameters().call()
    print('Optimization parameters received')
    L = np.array(L)
    L = L/10**6
    rho = rho/10**6
    if curtailPV==1: # sum(q_i)<0 <=> PVproduction > PVimported ==> curtailment of PV, no diminution of PV import
    	L[length:2*length] = np.minimum(L[length:2*length],0)
    	
    
    
    
    # re-define parameters for MATLAB
    rho_M =eng.double(rho)

    if len(Storage.optimizedProfile)==0:
        X_k_M = np.concatenate((np.zeros(length),np.abs(yPVforecast),np.zeros(2*length))) 
        X_k_M = eng.transpose(eng.cell2mat(X_k_M.tolist()))    
    else:
        X_k_M = Storage.optimizedProfile
        X_k_M = eng.transpose(eng.cell2mat(X_k_M.tolist()))
    #L_M = L[:,0]
    #L_M = L[0]
    L_M = eng.transpose(eng.cell2mat(L.tolist()))
    
    # optimization with MATLAB
    [X_storage, utility_value] = eng.util_func_B_4(beta1,beta2, parametersB,timestepM, startTimeM, endTimeM, price_PV_M, price_b_M, price_imp_M, rho_M, L_M, X_k_M, p_forecast_PV_total_M,nargout=2)
    X_storage = np.asarray(X_storage)
    Storage.optimizedProfile = X_storage[:,0]  
    profileToSend = tl.toListInt(Storage.optimizedProfile)
    tx_hash = contractADMM.functions.sendOptimizedProfile(profileToSend).transact()
    tx_receipt = web3.eth.waitForTransactionReceipt(tx_hash)
    gasSpend += tx_receipt.gasUsed
    print('Optimized profile send')
    

#    # optimization of profile
#    Storage.optimization(startTime, endTime, timeStep, price_imp, price_PV, price_B, rho, L, yPVforecast)     
#    profileToSend = tl.toListInt(Storage.optimizedProfile)
#    tx_hash = contractADMM.functions.sendOptimizedProfile(profileToSend).transact()
#    tx_receipt = web3.eth.waitForTransactionReceipt(tx_hash)
#    gasSpend += tx_receipt.gasUsed
#    print('Optimized profile send')

    # wait for global results
    l = 0
    while l!=1:
    #while l!=iterations:
        event_filter = contractADMM.events.EndIteration.createFilter(fromBlock='latest')
        #event_filter = contractADMM.events.EndIteration.createFilter(fromBlock=0)
        eventlist = event_filter.get_all_entries()
        l=len(eventlist)
    print('Global optimization ended')
    result = eventlist[0].args["_result"]
	
    iterations +=1		
    print('End of iteration')
    print(iterations)
    
    X_historic.append(Storage.optimizedProfile)
    L_historic.append(L)



#    # save data
    # SOC
    timestep = timeVector.freq
    timestep = timestep/np.timedelta64(1, 'h') # in hours 
    Storage.CalculateSOC(timestep)
    
    X_storage = Storage.optimizedProfile
    X_historic_storage = X_historic
    SOC_storage = Storage.SOC
    gasSpend_Storage = gasSpend
    gasSpend_without_initialization_Storage = gasSpend - gasSpend_initialization
    etherSpend_Storage = etherSpend - web3.eth.getBalance(web3.eth.defaultAccount)
    
    # save variables
    filename = 'Results_storage'
    with open(filename, 'wb') as f:
        pickle.dump(X_historic_storage,f)
        pickle.dump(SOC_storage,f)
        pickle.dump(gasSpend_Storage,f)
        pickle.dump(gasSpend_without_initialization_Storage,f)
        pickle.dump(etherSpend_Storage,f)
        pickle.dump(L_historic,f)


print('Final state')

web3.geth.miner.stop() 




# SOC
timestep = timeVector.freq
timestep = timestep/np.timedelta64(1, 'h') # in hours 
Storage.CalculateSOC(timestep)

X_storage = Storage.optimizedProfile
X_historic_storage = X_historic
SOC_storage = Storage.SOC
gasSpend_Storage = gasSpend
gasSpend_without_initialization_Storage = gasSpend - gasSpend_initialization
etherSpend_Storage = etherSpend - web3.eth.getBalance(web3.eth.defaultAccount)




# save variables
filename = 'Results_storage'
with open(filename, 'wb') as f:
    pickle.dump(X_historic_storage,f)
    pickle.dump(SOC_storage,f)
    pickle.dump(gasSpend_Storage,f)
    pickle.dump(gasSpend_without_initialization_Storage,f)
    pickle.dump(etherSpend_Storage,f)
    pickle.dump(L_historic,f)



# load the workspace:
#with open(filename, 'rb') as f:
 # X_historic =  pickle.load(f)
