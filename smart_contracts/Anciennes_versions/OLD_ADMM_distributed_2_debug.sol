pragma solidity ^0.5.1;

// As solidity only supports integers, we work with a factor 1 000 000.
// For example, for the real number 5.14, here in the smart contract we use 5.14*10**6


contract ADMM {
    uint N; // number of total users
    int Pmax; // maximal allower power
    uint nPV; // number of PV generators
    uint nPVforecast; // number of PV genertors who have send forecast
    uint numUsers;  // incremental number of users
	uint timeStep; // in secondes
	uint startTime; // in secondes (unix time)
	uint endTime;	// in secondes (unix time)
	uint nbSteps;  // number of steps for the vectors
	uint[] priceProfile; // dynamic size for the moment, to be changed
	int[] X_global; // sum of all the users' profile
	int[] Z; // global variable
	int[] Z_previous;
	int[] U; //global variable
    int[] L;
	int rho; // penalty parameter
	uint epsilon; // stop criteria
	int[] R;
	int[] S;
	uint normR;
	uint normS;
	int[] totalPVforecast; //
	uint finalResult; // 1 if final state reached, 0 else
	int tau; 
	uint mu;
	uint iterations;
	uint nbProfileSend;
	uint AAA; //////
	
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
//	event GlobalResults(address indexed _from, uint _result);
	event GlobalOptim(address indexed _from, uint _result);
	event EndIteration(address indexed _from, uint _result);
	event UserCreated(address indexed _from);
	
	
// ----------------------------------------------------------------------------------

	// constructor that registers the address of contract creator
	constructor(uint _N, int _Pmax, uint _epsilon) public {
	    Pmax = _Pmax;
	    N = _N;
	    epsilon = _epsilon;
		contractCreator = msg.sender;
		rho = 1*10**3;
		nPV = 0;
		nPVforecast = 0;
		finalResult = 0;
		normR = 0;
		normS = 0;
		mu = 10;
		tau =2; //10;
		iterations = 0;
		nbProfileSend = 0;
		AAA=0;////
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
	       Z.push(0);
	       Z_previous.push(0);
	       U.push(0);
	       L.push(0);
	       R.push(1*10**6); // inital value to enter into the optimization loop 
	       S.push(0);
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
	    
	    // send an event to inform that all users has updated his profile
	    nbProfileSend +=1;
	    if (nbProfileSend == N) {
	        emit SendProfile(msg.sender);
	    }
	    
	    
	    profileSend[msg.sender] =1; // prevent user to send twice a profile
	}
	
	
// ----------------------------------------------------------------------------------	

	function GlobalOptimization() public onlyOwner {
	    for (uint i = 0; i<nbSteps; i++) {
	        
	        Z_previous[i] = Z[i]; //
//	        if (iterations ==0) { // first iteration
//	            Z_previous[i] = X_global[i]/int(N);
//	        }
	        
	        // Z-update
	        if (X_global[i]>Pmax) {
	            Z[i] = Pmax/int(N) -3*10**(5);
	            //Z[i] = Pmax/N - 0.3;
	        }
	        else {
	          Z[i] = U[i] + X_global[i]/int(N); 
	        } 
	    }

		for (uint i = 0; i<nbSteps; i++) {
			// global variables update    
			U[i] = U[i] + X_global[i]/int(N)- Z[i];
			R[i] = X_global[i]/int(N) - Z[i];
			S[i] = rho*(Z[i] - Z_previous[i])/1000; //// attention, coeff différent pour rho
	    
	    // norm 2 of R and S
	    //normR += (abs(R[i]))**2;
	    //normS += (abs(S[i]))**2;
	    }
	    for (uint i = 0; i<nbSteps; i++) {
	        L[i] = X_global[i]/int(N) - Z[i] + U[i];
	    }
	    
	    normR = norm(R);
	    normS = norm(S);
	    
	    // stop criterion 
	    if (normR>epsilon**2 || normS>epsilon**2){
	        finalResult = 0;
	    } 
	    else {
	        finalResult = 1;
	    }
	    
	    // rho-update
	    if (normR>mu**2 * normS){
	        rho = tau*rho;
	        AAA = 1;
	    }
	    else if (normS>mu**2 * normR){
	        rho = rho/tau;
	        AAA = 2;
	    }
	    else{
	        AAA = 0;
	    }
	    
	    emit GlobalOptim(msg.sender,finalResult);
	}


//    function SendGlobalResults(int[] memory _L, int _rho, uint _result) public onlyOwner 	{
//		L = _L;
//		rho = _rho;
		// send event to inform users that optimization is ended
//        emit GlobalResults(msg.sender,_result);
//    }			
	
	
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
	    nbProfileSend = 0;
	    // send event to inform users that optimization is ended
	    emit EndIteration(msg.sender,finalResult);
    }

//    function InformUsersEnd() public onlyOwner {
//        // inform the users that global optimization ended
//         emit GlobalResults(msg.sender,finalResult);
//    }


	function getOptimizationParameters() public view returns (int[] memory, int){
	    return (L,rho);
	}

	

	
// ----------------------------------------------------------------------------------	
// TOOLS FUNCTION

    
    function abs(int _x) public pure returns (uint) {
        uint absValue;
        if (_x>0) {
            absValue = uint(_x);
        }
        else {
            absValue = uint(-_x);
        }
        return(absValue);
    }


    function norm(int[] memory _X) public pure returns (uint) {
        // returns norm(X)**2
        uint normX;
        normX=0;
        for (uint i = 0; i<_X.length; i++) {
            normX += (abs(_X[i]))**2;
        }
        return(normX);
    }



	function getX() public view returns (int[] memory){
	    return(X_global);
	}
	
	function getL() public view returns (int[] memory){
	    return(L);
	}
	
	function getRho() public view returns (int){
	    return(rho);
	}
	
	function getR() public view returns (int[] memory){
	    return(R);
	}
	
	function getS() public view returns (int[] memory){
	    return(S);
	}
	
	function getAAA() public view returns (uint){
	    return(AAA);
	}	

    function rhoUpdate(uint _normR, uint _normS, uint _mu, int _tau, int _rho) public pure returns (int){
        // rho-update
	    if (_normR>_mu**2 *_normS){
	        //_rho = _tau*_rho;
	        _rho = _tau*_rho;//
	    }
	    else if (_normS>_mu**2 * _normR){
	        _rho = _rho/_tau;
	    }
	    return(_rho);
    }

}