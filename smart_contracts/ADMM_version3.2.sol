pragma solidity ^0.6.0;

// As solidity only supports integers, we work with a factor 10**6
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
    uint[] priceProfile; // price of imported electricity (EPEX)
    uint coeffPricePV; // coeff price PV
    uint coeffPriceB; // coeff price battery
    int[] totalPVforecast; // sum of forecasted local PV production
    int[] sum_X_global; // sum of all the users' profile
    int[] meanX;
    int[] Z; // global variable
    int[] Z_previous;
    int[] U; //global variable
    int[] L;
	uint curtailPV;
    int rho; // penalty parameter
    uint epsilon; // stop criteria
    int[] R;
    int[] S;
    uint normR;
    uint normS;
    uint finalResult; // 1 if final state reached, 0 else
    int tau; 
    uint mu;
    uint iterations;
    uint nbProfileSend;
    
    address contractCreator;
    
    struct User {
	    uint userID;
	    uint class; // class (1:EV, 2: storage, 3: PV; 4: loads)
	    int[] optimProfile; // power profile of the user
	    //uint amount; // amount of money of the user
    }  
    
    address[] public userAccounts; // array that saves the address of users
    
    mapping(address => User) public users;
    mapping(address => uint) public registeredUsers; // (0: not registered, 1: already registered)
    mapping(address => uint) public profileSend; // (0: profile not send, 1: already send)
    mapping(address => uint) public PVforecastSend; // (0: profile not send, 1: already send)
    
    event SendProfile( address indexed _from);
    event SendForecast(address indexed _from); 
    event StartOptim(address indexed _from);
    //event GlobalResults(address indexed _from, uint _result);
    event GlobalOptim(address indexed _from, uint _result);
    event EndIteration(address indexed _from, uint _result);
    event UserCreated(address indexed _from);
	
	
// ----------------------------------------------------------------------------------

    // constructor that registers the address of contract creator
    constructor(uint _N, int _Pmax, uint _epsilon, int _rho) public {
    	Pmax = _Pmax;
    	N = _N;
	    epsilon = _epsilon;
	    contractCreator = msg.sender;
	    rho = _rho; //(1*10**(-6) *10**6= 1)
	    nPV = 0;
	    nPVforecast = 0;
	    finalResult = 0;
	    normR = 0;
	    normS = 0;
	    mu = 10;
	    tau = 5; //10;
	    iterations = 0;
	    nbProfileSend = 0;
		curtailPV = 2;
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
		for (uint i=0; i < 4*nbSteps; i++){ //4*nbSteps because 4 variables (p,q,lamda,r)
	       sum_X_global.push(0);
	       meanX.push(0);
	       Z.push(0);
	       Z_previous.push(0);
	       U.push(0);
	       L.push(0); 
	       R.push(1*10**6);  
	       S.push(0); 
	    }  
	    
		for (uint i=0; i < nbSteps; i++){
		    totalPVforecast.push(0);	
		}
		
	}
	
	
	function getTimeInterval() public view returns (uint, uint, uint){
		return (startTime,endTime,timeStep);
    }	
	
	
// ----------------------------------------------------------------------------------	
	function setPriceImp (uint[] memory _priceProfile)  public onlyOwner {
        // contractCreator can call this function to load the price vector for the desired period into the blockchain
		require(_priceProfile.length == nbSteps, "vector doesn't have good size");
	    for (uint i=0; i < nbSteps; i++){
	       priceProfile.push(_priceProfile[i]);
	   }
	}
	
	function setPriceLocal (uint  _coeffPricePV, uint _coeffPriceB) public onlyOwner {
	    // price_local = coeffPricePV* price_imp
	    coeffPricePV = _coeffPricePV;
	    coeffPriceB = _coeffPriceB;
        emit StartOptim(msg.sender); 
	   // once the agent has set time interval and price function, 
	   // this event indicates to the users to start the optimization
	    
	}
	
    function getPriceImp() public view returns (uint[] memory){
		return (priceProfile);
    }

    function getPriceLocal() public view returns (uint,uint){
        return (coeffPricePV, coeffPriceB);
    }
	
// ----------------------------------------------------------------------------------
		
	function sendOptimizedProfile(int[] memory _optimProfile) public {
	    require(_optimProfile.length == 4*nbSteps, "vector doesn't have good size");
	    require(profileSend[msg.sender]==0); // check that the user has not already send the profile 
	    User storage user = users[msg.sender];
	    user.optimProfile = _optimProfile;
	    for (uint i = 0; i<4*nbSteps; i++) {
	        sum_X_global[i] += user.optimProfile[i]; // update sum of profiles
	    }
	    
	    // send an event to inform that all users has updated his profile
	    nbProfileSend +=1;
	    if (nbProfileSend == N) {
	        emit SendProfile(msg.sender);
	    }
	    profileSend[msg.sender] =1; // prevent user to send twice a profile
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
// GLOBAL OPTIMIZATION

	function updateZ() public onlyOwner {
	    for (uint i = 0; i<4*nbSteps; i++) {
	        meanX[i] = sum_X_global[i]/int(N);
	        //if (iterations ==0) { // first iteration
	        //    Z[i] = sum_X_global[i]/int(N);
	        //}
	        Z_previous[i] = Z[i];
	    }
	    
	    // Z-update
	    for (uint i=0; i<nbSteps; i++) { // if global constraint violated
	        if (sum_X_global[2*nbSteps+i]>Pmax){ 
	            Z[i] = U[i] + meanX[i];   
	            Z[nbSteps+i] = 0;
	            Z[2*nbSteps+i] = Pmax/int(N) -3*10**(5);  //Z[i] = Pmax/N - 0.3;
	            Z[3*nbSteps+i] = 0;
	        }
	        else { // if global constraint respected
	            Z[i] = U[i] + meanX[i];
	            Z[nbSteps+i] = 0;
	            Z[2*nbSteps+i] = U[2*nbSteps+i] + meanX[2*nbSteps+i];
	            Z[3*nbSteps+i] = 0;
	        }
	    }
	}




    function updateURSrho() public onlyOwner {
        // update U, R, S and rho
	    for (uint i = 0; i<4*nbSteps; i++) {
			U[i] = U[i] + meanX[i]- Z[i];
			R[i] = meanX[i] - Z[i];
			S[i] = rho*(Z[i] - Z_previous[i])/10**(6);//
			
		}
		normR = norm(R);
		normS = norm(S);
		
		
		// update rho
	    if (normR>mu**2 * normS){
	        rho = tau*rho;
	    }
	    else if (normS>mu**2 * normR){
	        rho = rho/tau;
	    }
	    
        // stop criterion 
	    if (normR>epsilon**2 || normS>epsilon**2){
	        finalResult = 0;
	    } 
	    else {
	        finalResult = 1;
	    }
    }
    
    
    
    function updateL() public onlyOwner {
        // update L (L_PV and L_2)
        for (uint i = 0; i<4*nbSteps; i++) {
            L[i] = meanX[i] - Z[i] + U[i];
		}
            
        
		if(getArraySum(sum_X_global, nbSteps, 2*nbSteps)>0) { //sum(q_i)>0 --> importsPV > productionPV ==> no PV curtailment, diminution of PV imports from consumers
			curtailPV = 0;
		}
        else{//sum(q_i)<0 --> importsPV < productionPV ==> no diminution of PV import but PV curtailment 
			curtailPV = 1;
        }
        
        emit GlobalOptim(msg.sender,finalResult);
    }
    
    
    
	
	
	function reset() public onlyOwner{
        // reset variables at each iterations steps
        // (mappings, rho, L, X, forecast...)
        for (uint i = 0; i<N; i++) {
            
            address _addr = userAccounts[i];
	        profileSend[_addr] = 0; 
        }
	    for (uint i = 0; i<4*nbSteps; i++) {
	        sum_X_global[i] = 0; 
	    }
	    nbProfileSend = 0;
	    iterations +=1;
	    // send event to inform users that optimization is ended
	    emit EndIteration(msg.sender,finalResult);
    }



	function getOptimizationParameters() public view returns (int[] memory, int, uint){
	    return (L,rho,curtailPV);
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


    function getArraySum(int[]  memory _array, uint _start, uint _end)  public pure returns (int sum_) 
    {
        sum_ = 0;
        for (uint i = _start; i < _end; i++) {
            sum_ += _array[i];
        }
    }



	function get_sumX() public view returns (int[] memory){
	    return(sum_X_global);
	}
	
	function getL() public view returns (int[] memory){
	    return(L);
	}
	
	function getZ() public view returns (int[] memory){
	    return(Z);
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
	
	function getnormR() public view returns (uint){
	    return(normR);
	}
	
	function getnormS() public view returns (uint){
	    return(normS);
	}
	
	
}