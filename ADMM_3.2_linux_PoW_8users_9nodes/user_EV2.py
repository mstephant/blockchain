# -*- coding: utf-8 -*-

"""
Last modification: 15/09/2020
@author: matthieu.stephant

script for EV profile optimization 
interaction with blockchain node
PoW
smart contract: ADMM_version3.2.sol
"""


import json
import classes_definition as cl
import numpy as np
import pandas as pd
import tools as tl 
import pickle
#from web3.middleware import geth_poa_middleware
import time

#from web3.auto import web3 # IPC automatic connection
from web3 import Web3
node_url = "http://localhost:8001"
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
EV2 = cl.EV()

# gas counter initialization
gasSpend = 0
etherSpend = web3.eth.getBalance(web3.eth.defaultAccount)

# create a User struct
tx_hash = contractADMM.functions.newUser(EV2.userType).transact()
print('Call newUser function')

tx_receipt = web3.eth.waitForTransactionReceipt(tx_hash)
print('newUser transaction mined')
gasSpend += tx_receipt.gasUsed

# preferences
EV2.setPreferences(1/3,1/3,1/3) # (costs, exchanges with community, comfort)
    
# parameters
PmaxEV = 22
EmaxEV = 41
etaEV = 80/100
SOCmin = 20/100
SOCmax = 100/100
SOCinit = 45/100
SOCrequired = 80/100
arrivalTimeEV2 = pd.to_datetime('2019-06-10 10:30')
departureTimeEV2 = pd.to_datetime('2019-06-10 19:00')

EV2.setChargingTime(arrivalTimeEV2,departureTimeEV2)
EV2.setParameters(PmaxEV,EmaxEV,etaEV,SOCmin,SOCmax,SOCinit)
EV2.setSOCrequired(SOCrequired)


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
L_historic = []

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
    	

    # optimization of profile
    EV2.optimization(startTime, endTime, timeStep, price_imp, price_PV, price_B, rho, L, yPVforecast) 
    profileToSend = tl.toListInt(EV2.optimizedProfile)
    tx_hash = contractADMM.functions.sendOptimizedProfile(profileToSend).transact()
    tx_receipt = web3.eth.waitForTransactionReceipt(tx_hash)
    gasSpend += tx_receipt.gasUsed
    print('Optimized profile send')

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
    
    X_historic.append(EV2.optimizedProfile)
    
    L_historic.append(L)
    
#    # save data
    # SOC
    timestep = timeVector.freq
    timestep = timestep/np.timedelta64(1, 'h') # in hours 
    EV2.CalculateSOC(timestep)
    
    X_EV2 = EV2.optimizedProfile
    X_historic_EV2 = X_historic
    SOC_EV2 = EV2.SOC
    gasSpend_EV2 = gasSpend
    gasSpend_without_initialization_EV2 = gasSpend - gasSpend_initialization
    etherSpend_EV2 = etherSpend - web3.eth.getBalance(web3.eth.defaultAccount)
    
    # save variables
    filename = 'Results_EV2'
    with open(filename, 'wb') as f:
        pickle.dump(X_historic_EV2,f)
        pickle.dump(SOC_EV2,f)
        pickle.dump(gasSpend_EV2,f)
        pickle.dump(gasSpend_without_initialization_EV2,f)
        pickle.dump(etherSpend_EV2,f)
        pickle.dump(L_historic,f)

    
print('Final state')
web3.geth.miner.stop()



# SOC
timestep = timeVector.freq
timestep = timestep/np.timedelta64(1, 'h') # in hours 
EV2.CalculateSOC(timestep)

X_EV2 = EV2.optimizedProfile
X_historic_EV2 = X_historic
SOC_EV2 = EV2.SOC
gasSpend_EV2 = gasSpend
gasSpend_without_initialization_EV2 = gasSpend - gasSpend_initialization
etherSpend_EV2 = etherSpend - web3.eth.getBalance(web3.eth.defaultAccount)





# save variables
filename = 'Results_EV2'
with open(filename, 'wb') as f:
    pickle.dump(X_historic_EV2,f)
    pickle.dump(SOC_EV2,f)
    pickle.dump(gasSpend_EV2,f)
    pickle.dump(gasSpend_without_initialization_EV2,f)
    pickle.dump(etherSpend_EV2,f)
    pickle.dump(L_historic,f)
    


# load the workspace:
#with open(filename, 'rb') as f:
 # X_historic =  pickle.load(f)
