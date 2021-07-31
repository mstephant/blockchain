pragma solidity ^0.5.1;

// As solidity only supports integers, we work with a factor 1 000 000.
// For example, for the real number 5.14, here in the smart contract we use 5.14*10**6


contract ADMM {
    uint N; // number of total users
    uint nPV; // number of PV genertors
    uint nPVforecast; // number of PV genertors who have send forecast
    uint numUsers;  // incremental number of users
	uint timeStep; // in secondes
	uint startTime; // in secondes (unix time)
	uint endTime;	// in secondes (unix time)
	uint nbSteps;  // number of steps for the vectors
	uint[] priceProfile; // dynamic size for the moment, to be changed
	int[] X_global; // sum of all the users' profile
    int[] L;
	uint rho; // penalty parameter
	int[] totalPVforecast; //
	
	address contractCreator;

	struct User {
	    uint userID;
		uint class; // class (1:EV, 2: storage, 3: PV; 4: loads)
		int[] optimProfile; // power profile of the user
		//uint amount; // amount of money of the user
		//int[] forecastProfile; // dynamic size for the moment, to be changed // EST-CE UTILE ? (seul l'utilisateur doit connaître localement son profil prévu)
	}
	
	address[] public userAccounts; // array that saves the address of users

	
	mapping(address => User) public users;
	mapping(address => uint) public registeredUsers; // (0: not registered, 1: already registered)
	mapping(address => uint) public profileSend; // (0: profile not send, 1: already send)
	mapping(address => uint) public PVforecastSend; // (0: profile not send, 1: already send)
	
		
	event SendProfile( address indexed _from);
	event SendForecast(address indexed _from); 
	event StartOptim(address indexed _from);
	event GlobalResults(address indexed _from, uint _result);
	event UserCreated(address indexed _from);
	    
// ----------------------------------------------------------------------------------

	// constructor that registers the address of contract creator
	constructor(uint _N) public {
	    N = _N;
		contractCreator = msg.sender;
		rho = 1;
		nPV = 0;
		nPVforecast = 0;
	}
	
	// modifier that checks if the person calling the function is the owner
	// (if the account calling the function is not the owner, then Solidity will throw an error)
	modifier onlyOwner() {
		require(msg.sender == contractCreator);
    	_;
	}

// ----------------------------------------------------------------------------------
    
    function newUser(uint _class) public returns (uint _userID){
	// user calls this function to declare its class
    // each user gets a different _userID	
        require(registeredUsers[msg.sender]==0); // check that the user is not already registered
	    _userID = numUsers++;
	    int[] memory _zero;
	    users[msg.sender] = User(_userID,_class, _zero);
	    registeredUsers[msg.sender] =1; // prevent user to create twice a user struct
	    userAccounts.push(msg.sender);
	    if (_class ==3) { // count the number of PV generators
	        nPV +=1;
	    }
	    else {}
	    
	    emit UserCreated(msg.sender); // inform that a User has been created
	}		
    
  
// ----------------------------------------------------------------------------------
	function setTimeInterval(uint _startTime, uint _endTime, uint _timeStep) public onlyOwner {
	    // chose carefully startTime, endTime and timeStep
		require(_startTime < _endTime, "startTime shouble be less than endTime");
		require(_timeStep < _endTime - _startTime, "impossible timestep");
		nbSteps = (_endTime - _startTime)/_timeStep+1; // length of vector time
		// convertir en datetime ? 
		startTime = _startTime;
		endTime = _endTime;
		timeStep = _timeStep;
		
		// create zero vector with good length
		for (uint i=0; i < nbSteps; i++){
	       X_global.push(0);
	       L.push(0);
	       totalPVforecast.push(0);
	   }   
	}
	
	
	function getTimeInterval() public view returns (uint, uint, uint){
		return (startTime,endTime,timeStep);
    }	
	
	
// ----------------------------------------------------------------------------------	
	function setPrice (uint[] memory _priceProfile)  public onlyOwner {
        // contractCreator can call this function to load the price vector for the desired period into the blockchain
		require(_priceProfile.length == nbSteps, "vector doesn't have good size");
	    for (uint i=0; i < nbSteps; i++){
	       priceProfile.push(_priceProfile[i]);
	   }
	   emit StartOptim(msg.sender); 
	   // once the agent has set time interval and price function, 
	   // this event indicates to the users to start the optimization
	}
	
	
    function getPrice() public view returns (uint[] memory){
		return (priceProfile);
    }

	
// ----------------------------------------------------------------------------------
	function sendPVforecast(int[] memory _PVforecast) public {
	// PV user send its forecast 
		require(_PVforecast.length == nbSteps, "vector doesn't have good size");
		User storage user = users[msg.sender];
	    require(user.class == 3); // only PV generator can send profile
	    require(PVforecastSend[msg.sender]==0); // check that the user has not already send the profile 
	    for (uint i = 0; i<nbSteps; i++) {
	        totalPVforecast[i] += _PVforecast[i];
	    }
	    nPVforecast +=1; // increment number of PV that hav send forecast
	    PVforecastSend[msg.sender] =1; // prevent user to send twice a profile
	    // send an event to inform that all PV generators have send their forecast
	    if (nPVforecast == nPV) {
	       emit SendForecast(msg.sender); 
	    }
	    else {}
	}
	
	
	function getPVforecast() public view returns (int[] memory){
	    return (totalPVforecast);
	}

// ----------------------------------------------------------------------------------	
		
	function sendOptimizedProfile(int[] memory _optimProfile) public {
	    require(_optimProfile.length == nbSteps, "vector doesn't have good size");
	    require(profileSend[msg.sender]==0); // check that the user has not already send the profile 
	    User storage user = users[msg.sender];
	    user.optimProfile = _optimProfile;
	    for (uint i = 0; i<nbSteps; i++) {
//	        if (user.class == 3) {
//	            user.optimProfile[i] = -_optimProfile[i]; // negative values for PV generators
//	        }
	        X_global[i] += user.optimProfile[i]; // update sum of profiles
	    }
	    // send an event to inform that user has updated his profile
	    emit SendProfile(msg.sender);
	    profileSend[msg.sender] =1; // prevent user to send twice a profile
	}
	
	
// ----------------------------------------------------------------------------------	

	function getGlobalResults() public onlyOwner view returns(int[] memory) {
	    return(X_global);   
	}


    function SendGlobalResults(int[] memory _L, uint _rho, uint _result) public onlyOwner 	{
		L = _L;
		rho = _rho;
		// send event to inform users that optimization is ended
        emit GlobalResults(msg.sender,_result);
    }			
	
	function getOptimizationParameters() public view returns (int[] memory, uint){
	    return (L,rho);
	}

	
    function reset() public onlyOwner{
        // reset variables at each iterations steps
        // (mappings, rho, L, X, forecast...)
        for (uint i = 0; i<N; i++) {
            address _addr = userAccounts[i];
	        profileSend[_addr] = 0; 
        }
	    for (uint i = 0; i<nbSteps; i++) {
	        X_global[i] = 0;
	    }
    }
	
	
// ----------------------------------------------------------------------------------	
// FONCTIONS DE TEST

	function getX() public view returns (int[] memory){
	    return(X_global);
	}
	
	
	function getAllAddresses() public onlyOwner view returns (address[]  memory){
	    return(userAccounts);
	}
	
	function getNPV() public view returns (uint) {
	    return(nPV);
	}




//	function moneyTransfer {
//	}


}