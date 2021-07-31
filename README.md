# Local blockchain for distributed optimization on a local energy community

Test of energy consumption of a local Ethereum blockchain for distributed optimization of energy on a local energy community, containing diverse actors (electric vehicles, PV generators, tertiary loads, storage system).

The smart contract implements a distributed optmization algorithm that enables the actors to find a consensus on the share of the local PV generation. 

The python files are used to implement the behavior of the actor and to automate the exchange with blockchain nodes. In our model, each user holds one node of the blockchain.

We test the implementation with PoA and PoW. 

Tested on March 2021 with Geth under Linux distribution


Steps to follow:
- Start a bootnode:
- 
  bootnode --genkey=boot.key
  
  bootnode --nodekey=boot.key
  
  
- Set up one node for each user of the community, plus one additonnal, with the following commands (use puppeth to generate the genesis file):
  
  geth --datadir node0 init genesis.json
  
  geth --datadir "node0" --identity "Node0" --rpc --rpcport "8000" --rpccorsdomain "*" --port "30303" --rpcapi "eth,net,web3,miner,personal,debug " --networkid \<networkid value\> --syncmode "full" --allow-insecure-unlock --bootnodes \<enode value\>  console

  
- Start the python scripts starting by 'user_...'. Modify if needed the address of the node in order to make sure that the connection is successful. Each python script contains some physical parameters of the actor it represents, that can be modified. 


