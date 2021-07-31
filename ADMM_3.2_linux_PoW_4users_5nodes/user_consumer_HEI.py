# -*- coding: utf-8 -*-

"""
Last modification: 23/09/2020
@author: matthieu.stephant

script for load
interaction with blockchain node
PoA with 1 node only
smart contract:ADMM_version3.2.sol
	-> calculation of newL is done in the python code and not in the smart contract

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
node_url = "http://localhost:8004"
web3 = Web3(Web3.HTTPProvider(node_url,request_kwargs={'timeout':1000}))

# set pre-funded account as sender
web3.eth.defaultAccount = web3.eth.accounts[0]
web3.geth.personal.unlockAccount(web3.eth.defaultAccount,"0000",0) # unlock account permanently

# mining of all sealer nodes necessary
web3.geth.miner.start()

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
Building = cl.Consumer()

# gas counter initialization
gasSpend = 0
etherSpend = web3.eth.getBalance(web3.eth.defaultAccount)

# create a User struct
tx_hash = contractADMM.functions.newUser(Building.userType).transact()
tx_receipt = web3.eth.waitForTransactionReceipt(tx_hash)
gasSpend += tx_receipt.gasUsed

# preferences
Building.setPreferences(1/3,1/3,1/3)  #(costs, local exchanges,comfort)    
# parameters
flexibility = 25/100
Building.setFlexibility(flexibility)


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

if startTime.month == 6:
       month = 'juin'
elif startTime.month == 12:
       month = 'decembre'
       
       


# get price vector 
price_imp = contractADMM.functions.getPriceImp().call() 
price_imp = np.array(price_imp)
price_imp = pd.Series(price_imp/10**6,timeVector)

[coeff_pricePV, coeff_priceB] = contractADMM.functions.getPriceLocal().call() 
coeff_pricePV = coeff_pricePV/10**6
price_PV = pd.Series(data=coeff_pricePV*np.ones(length),index=timeVector)
coeff_priceB = coeff_priceB/10**6
price_B = pd.Series(data=coeff_priceB*np.ones(length),index=timeVector)

# get forecast_consumer_HEI profile
fnamePV ='../Data_' + month + '19/conso_HEI_' + month + '19.csv'
forecast_consumer_HEI =  pd.read_csv(fnamePV,sep=';',names = ['timeVector','power'],index_col=0,squeeze=True,parse_dates=True)
forecast_consumer_HEI = forecast_consumer_HEI.asfreq(freq =str(timeStep)+'T') # retime
forecast_consumer_HEI = forecast_consumer_HEI.interpolate()
forecast_consumer_HEI = forecast_consumer_HEI[startTime:endTime]
Building.setForecast(forecast_consumer_HEI)
forecast_consumer_HEI = tl.toListInt(forecast_consumer_HEI)

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
    Building.optimization(startTime, endTime, timeStep, price_PV, price_imp, price_B, rho, L, yPVforecast)
    profileToSend = tl.toListInt(Building.optimizedProfile)
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
    
    X_historic.append(Building.optimizedProfile)
    L_historic.append(L)
    
    
    
#    #save date
    X_consumer_HEI = Building.optimizedProfile
    X_historic_consumer_HEI = X_historic
    gasSpend_consumer_HEI = gasSpend
    gasSpend_without_initialization_consumer_HEI = gasSpend - gasSpend_initialization
    etherSpend_consumer_HEI = etherSpend - web3.eth.getBalance(web3.eth.defaultAccount)
    #forecast_consumer_HEI = pd.Series(np.array(forecast_consumer_HEI)/10**6,timeVector)
    forecast_consumer_HEI = np.array(forecast_consumer_HEI)/10**6
    
    
    # save variable X_historic
    filename = 'Results_consumer_HEI'
    with open(filename, 'wb') as f:
        pickle.dump(X_historic_consumer_HEI,f)
        pickle.dump(forecast_consumer_HEI,f)
        pickle.dump(gasSpend_consumer_HEI,f)
        pickle.dump(gasSpend_without_initialization_consumer_HEI,f)
        pickle.dump(etherSpend_consumer_HEI,f)
        pickle.dump(L_historic,f)
    


print('Final state')
web3.geth.miner.stop()

X_consumer_HEI = Building.optimizedProfile
X_historic_consumer_HEI = X_historic
gasSpend_consumer_HEI = gasSpend
gasSpend_without_initialization_consumer_HEI = gasSpend - gasSpend_initialization
etherSpend_consumer_HEI = etherSpend - web3.eth.getBalance(web3.eth.defaultAccount)
#forecast_consumer_HEI = pd.Series(np.array(forecast_consumer_HEI)/10**6,timeVector)
forecast_consumer_HEI = np.array(forecast_consumer_HEI)/10**6


# save variable X_historic
filename = 'Results_consumer_HEI'
with open(filename, 'wb') as f:
    pickle.dump(X_historic_consumer_HEI,f)
    pickle.dump(forecast_consumer_HEI,f)
    pickle.dump(gasSpend_consumer_HEI,f)
    pickle.dump(gasSpend_without_initialization_consumer_HEI,f)
    pickle.dump(etherSpend_consumer_HEI,f)
    pickle.dump(L_historic,f)
    
    



# load the workspace:
#with open(filename, 'rb') as f:
 # X_historic =  pickle.load(f)
