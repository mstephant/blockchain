# -*- coding: utf-8 -*-

"""
Last modification: 23/09/2020
@author: matthieu.stephant

script for agent client
interaction with blockchain node
PoW
smart contract: ADMM_version3.2.sol
	-> calculation of newL is done in the python code and not in the smart contract

"""

import json
import pandas as pd
import tools as tl
import numpy as np
import matplotlib.pyplot as plt 
import pickle
#from web3.middleware import geth_poa_middleware

#from web3.auto import web3 # IPC automatic connection
from web3 import Web3
node_url = "http://localhost:8000"
web3 = Web3(Web3.HTTPProvider(node_url,request_kwargs={'timeout':5000}))

# PoA compatibility
#web3.middleware_onion.inject(geth_poa_middleware, layer=0)

# set pre-funded account as sender
web3.eth.defaultAccount = web3.eth.accounts[0]
web3.geth.personal.unlockAccount(web3.eth.defaultAccount,"0000",0) # unlock account permanently

# mining of all sealers nodes necessary
web3.geth.miner.start()

# gas counter initialization
gasSpend = 0
etherSpend = web3.eth.getBalance(web3.eth.defaultAccount)

# parameters
N = 8 # number of players
Pmax = 350 # max allowed power (mW)
epsilon = 1 # stop criterion
rho = 1*10**(-4)

LCOE_PV = 0.8 # prix_local = coeff_prix * prix_imp
LCOS_B = 0.8;


### DEPLOY CONTRACT
#import contract details
with open("../../smart_contracts/ADMM_version3.2.json",'r') as f:
    abi = json.load(f)    

with open("../../smart_contracts/ADMM_version3.2.bin",'r') as f:
    bytecode = f.read()

# set pre-funded account as sender
web3.eth.defaultAccount = web3.eth.accounts[0]
#web3.geth.personal.unlockAccount(web3.eth.defaultAccount,"0000",0) # unlock account permanently


### CONTRACT DEPLOYMENT
# Instantiate and deploy contract
ADMM = web3.eth.contract(abi=abi, bytecode=bytecode)
tx_hash = ADMM.constructor(N,int(Pmax*10**6),int(epsilon*10**6),int(rho*10**6)).transact() 

# Wait for the transaction to be mined, and get the transaction receipt
tx_receipt = web3.eth.waitForTransactionReceipt(tx_hash)
gasSpend += tx_receipt.gasUsed
print(tx_receipt.contractAddress)

# write address in a seperate file
f = open("address.txt", "w")
f.write(tx_receipt.contractAddress)
f.close()

# Create the contract instance with the newly-deployed address
contractADMM = web3.eth.contract(
    address=tx_receipt.contractAddress,
    abi=abi,
)
print('Contract deployed')


### ENERGY COMMUNITY MANAGER
# wait that all users create struct User in smart contract
l = 0
while l!=N:
    event_filter = contractADMM.events.UserCreated.createFilter(fromBlock=0)
    eventlist = event_filter.get_all_entries()
    l=len(eventlist)
print('All users have created struct in smart contract')



# specify time interval
startTime = pd.to_datetime('2019-06-10 00:00')  # one day in june or december
endTime = pd.to_datetime('2019-06-11 00:00')    # one day in june or december
timeStep = 15 # minutes
timeVector = pd.date_range(startTime, endTime,freq = str(timeStep)+'T')
length = len(timeVector)

if startTime.month == 6:
       month = 'juin'
elif startTime.month == 12:
       month = 'decembre'


# get price function
fname ='../Data_' + month + '19/price_' + month + '19.csv'
price_imp = pd.read_csv(fname,sep=';',names = ['timeVector','price'],index_col=0,squeeze=True,parse_dates=True)
price_imp = price_imp.asfreq(freq =str(timeStep)+'T', method = 'ffill')
price_imp = 1/1000* price_imp[startTime:endTime]  # from €/Mwh to €/kWh
price_imp = tl.toListInt(price_imp)

price_PV = LCOE_PV*np.ones(length)
price_B = LCOS_B*np.ones(length)

# call setTimeInterval
startTime = int(startTime.timestamp())
endTime = int(endTime.timestamp())
timeStep = timeStep*60
tx_hash = contractADMM.functions.setTimeInterval(startTime,endTime,timeStep).transact()
tx_receipt = web3.eth.waitForTransactionReceipt(tx_hash)
gasSpend += tx_receipt.gasUsed
print('Set time')

# call setPriceImp and setPriceLocal
tx_hash = contractADMM.functions.setPriceImp(price_imp).transact()
tx_receipt = web3.eth.waitForTransactionReceipt(tx_hash)
gasSpend += tx_receipt.gasUsed

tx_hash = contractADMM.functions.setPriceLocal(int(LCOE_PV*10**6),int(LCOS_B*10**6)).transact() 
tx_receipt = web3.eth.waitForTransactionReceipt(tx_hash)
gasSpend += tx_receipt.gasUsed
print('Set price')


# intialization variables
iterations = 0
result = 0
number_last_block = web3.eth.getBlock('latest').number


# variables to follow convergence
L_historic=[]
X_historic=[]
R_historic = []
S_historic = []
rho_historic = []
Z_historic = []

gasSpend_initialization = gasSpend

while result == 0:
    
    # wait for users to update their profile
    l = 0
    print('Waiting for users profiles')
    while l!=1:
        new_latest = web3.eth.getBlock('latest').number
        if new_latest > number_last_block:
            number_last_block+=1
            event_filter = contractADMM.events.SendProfile.createFilter(fromBlock=number_last_block)
            eventlist = event_filter.get_all_entries()
            l=len(eventlist) 
    print('Received profiles from all users')
    
    
    # process global optimization
    print("Starting global optimization")
    print("Z update")
    tx_hash = contractADMM.functions.updateZ().transact()
    tx_receipt = web3.eth.waitForTransactionReceipt(tx_hash)
    gasSpend += tx_receipt.gasUsed
    
    number_last_block = web3.eth.getBlock('latest').number
    
    print("U, R, S, L and rho update")
    tx_hash = contractADMM.functions.updateURSrho().transact()
    tx_receipt = web3.eth.waitForTransactionReceipt(tx_hash)
    gasSpend += tx_receipt.gasUsed
    
    tx_hash = contractADMM.functions.updateL().transact()
    tx_receipt = web3.eth.waitForTransactionReceipt(tx_hash)
    gasSpend += tx_receipt.gasUsed
    
    
    # wait for global results
    l = 0
#    while l!=iterations:
#        event_filter = contractADMM.events.GlobalOptim.createFilter(fromBlock=0)
    while l!=1:
        event_filter = contractADMM.events.GlobalOptim.createFilter(fromBlock='latest')
        eventlist = event_filter.get_all_entries()
        l=len(eventlist)
    print('Global optimization ended')
    
    
    
    ### check values
    Z = np.array(contractADMM.functions.getZ().call())
    Z = np.array(Z/10**6)
    Z_historic.append(Z)
    X = np.array(contractADMM.functions.get_sumX().call())
    X = np.array(X/10**6)
    X_historic.append(X)
    L = np.array(contractADMM.functions.getL().call())
    L = np.array(L/10**6)
    L_historic.append(L)
    R = np.array(contractADMM.functions.getR().call())
    R = np.array(R/10**6)
    R_historic.append(R)
    S = np.array(contractADMM.functions.getS().call())
    S = np.array(S/10**6)
    S_historic.append(S)
    rho = np.array(contractADMM.functions.getRho().call())
    rho_historic.append(np.array(rho*10**(-6)))
    normR = np.array(contractADMM.functions.getnormR().call())
    normS = np.array(contractADMM.functions.getnormS().call())

    
    # reset lists in smart contract
    tx_hash = contractADMM.functions.reset().transact() # reset data in the smart contract before new iteration
    tx_receipt = web3.eth.waitForTransactionReceipt(tx_hash)
    gasSpend += tx_receipt.gasUsed
    
    
    result = eventlist[0].args["_result"] 
    print('Rho value:')
    print(rho*10**(-6))
    print('normR:')
    print(np.sqrt(normR)*10**(-6))
    print('normS:')
    print(np.sqrt(normS)*10**(-6))
    iterations +=1
    print('End of iterations, reset values')
    print(iterations)





    gasSpend_agent = gasSpend
    gasSpend_without_initialization_agent = gasSpend - gasSpend_initialization
    etherSpend_agent = etherSpend - web3.eth.getBalance(web3.eth.defaultAccount)
    
    # save variables
    filename = 'Results_Agent'
    with open(filename, 'wb') as f:
        pickle.dump(timeVector,f)
        pickle.dump(price_imp,f)
        pickle.dump(LCOE_PV,f)
        pickle.dump(LCOS_B,f)
        pickle.dump(X_historic,f)
        pickle.dump(L_historic,f)
        pickle.dump(R_historic,f)
        pickle.dump(S_historic,f)
        pickle.dump(rho_historic,f)
        pickle.dump(iterations,f)
        pickle.dump(Pmax,f)
        pickle.dump(gasSpend_agent,f)
        pickle.dump(gasSpend_without_initialization_agent,f)
        pickle.dump(etherSpend_agent,f)


print('Final state')

web3.geth.miner.stop()

gasSpend_agent = gasSpend
gasSpend_without_initialization_agent = gasSpend - gasSpend_initialization
etherSpend_agent = etherSpend - web3.eth.getBalance(web3.eth.defaultAccount)

# save variables
filename = 'Results_Agent'
with open(filename, 'wb') as f:
    pickle.dump(timeVector,f)
    pickle.dump(price_imp,f)
    pickle.dump(LCOE_PV,f)
    pickle.dump(LCOS_B,f)
    pickle.dump(X_historic,f)
    pickle.dump(L_historic,f)
    pickle.dump(R_historic,f)
    pickle.dump(S_historic,f)
    pickle.dump(rho_historic,f)
    pickle.dump(iterations,f)
    pickle.dump(Pmax,f)
    pickle.dump(gasSpend_agent,f)
    pickle.dump(gasSpend_without_initialization_agent,f)
    pickle.dump(etherSpend_agent,f)
    
# load the workspace:
#with open(filename, 'rb') as f:
 # X_historic =  pickle.load(f)
