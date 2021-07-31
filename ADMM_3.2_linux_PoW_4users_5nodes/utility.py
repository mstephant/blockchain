# -*- coding: utf-8 -*-

"""
Last modification: 15/09/2020
@author: matthieu.stephant

Optimization functions for all classes (EV, storage, PV generator, load)
"""

import numpy as np
import pandas as pd
import scipy.optimize as opt



###############################################################################

def optimEV(alpha, arrivalTime, departureTime, parametersEV,SOCrequired,\
                   startTime, endTime, timestep, price_imp, price_PV, price_B,\
                   rho, L, yK, yPVforecast):
    
    """ optimization of the EV utility function with penalty parameter rho        
        
        U_EV(x) = -alpha_1*(pi_PV*q+pi_imp*lambda+pi_b*r)
                 - alpha_2*(p_forecast_PV_total - q)**2
                 + alpha_3*ln(1-p)     
        
        function to minimize: -sum(U_EV(x)) + (rho/2)*||x-yK+L||**2
        

        with variables:
            p: net production (<0 if consumption)
            q: energy imported from the energy community (>0 if imported, <0 if exported)
            lambda: energy imported from external grid (>0)
            r: energy imported from the battery (>0 if imported, <0 if exported)
        
        
        
        Inputs: 
            - alpha = [alpha1,alpha2,alpha3] (costs, exchanges with community, comfort)
            - parametersEV = [PmaxEV,EmaxEV,etaEV,SOCminEV,SOCmaxEV,SOCinitEV]
            - SOCrequired: float
            - arrivalTime, departureTime: pd.to_datetime()
            - startTime, endTime: pd.to_datetime()
            - timestep: minutes
            - price_imp, price_PV, price_B, yPVforecast: pd.Series
            - yK, L: np.array
            - rho: floats
            - yPVforecast: pd.Series, total PV forecast
            
        Outputs:
            - yEV: np.array - variables (p,q,lambda) that maximizes utitility               
    """
    
    PmaxEV = parametersEV[0]
    EmaxEV = parametersEV[1]
    etaEV = parametersEV[2]
    SOCmin = parametersEV[3]
    SOCmax = parametersEV[4]
    SOCinit = parametersEV[5]
    
    timeVector = pd.date_range(startTime, endTime,freq = str(timestep)+'T')
    length = len(timeVector) # vector length
    
    # selection on the correct time interval
#    price_imp = price_imp[startTime:endTime] 
#    price_PV = price_PV[startTime:endTime]
#    price_B = price_B[startTime:endTime]

    price_imp = np.array(price_imp[startTime:endTime])
    price_PV = np.array(price_PV[startTime:endTime])
    price_B = np.array(price_B[startTime:endTime])
    

   
    yPVforecast = yPVforecast[startTime:endTime]
    yPVforecast = np.array(np.abs(yPVforecast))
    yPVforecast = yPVforecast+0.000001 # avoid 0 as it creates problem that I don't understand
    
    
    # intial step
    if len(yK) !=0: 
        pass
    else:
        #yK=np.zeros(4*length) # initial step
        yK = np.concatenate((np.zeros(length),yPVforecast,np.zeros(length),np.zeros(length)))
        


    # utility function to optimize
    def utilityEV(x):
        """
        U_EV(x) = -alpha_1*(pi_PV*q+pi_imp*lambda+pi_b*r)
                 - alpha_2*(p_forecast_PV_total - q)**2
                 + alpha_3*ln(1-p)     
        
        function to minimize: -sum(U_EV(x)) + (rho/2)*||x-yK+L||**2
        
        p = x[0:length]
        q = x[length:2*length]
        lambda = x[2*length:3*length]
        r = x[3*length:4*length]
        """
        utility = - alpha[0]*(price_PV*x[length:2*length]+price_imp*x[2*length:3*length]+price_B*x[3*length:4*length])/((np.maximum.reduce([price_imp, price_PV, price_B]))*PmaxEV) \
                  - alpha[1]*(yPVforecast-x[length:2*length])**2 /np.max(yPVforecast-PmaxEV)**2 \
                  + alpha[2]*np.log(1-x[0:length])/np.log(1+PmaxEV)
 

        sumUtility = -np.sum(utility) + (rho/2)*np.linalg.norm(x-yK+L)**2   
   
        return(sumUtility)
      
    
    # bounds: 0 <= |p(t)| <= Pmax
    LB_p = -PmaxEV*np.ones(length)
    UB_p = 1.0*np.zeros(length)
    
    # 0<q< yPVforecast (q>0 because EV can not export energy)
    LB_q = 1.0*np.zeros(length)
    UB_q = yPVforecast
    
    # lambda >0 and r >0
    LB_lambda_r = 1.0*np.zeros(2*length)
    UB_lambda_r = np.inf*np.ones(2*length)
     
    LB = np.concatenate((LB_p,LB_q,LB_lambda_r))
    UB = np.concatenate((UB_p,UB_q,UB_lambda_r))
    bnds = opt.Bounds(LB,UB)

    # equality constraints: 
    eq_constraints = []
    
        # p(t)=0 when EV not connected to the station
    a = pd.Series(np.ones(length),timeVector)
    a[arrivalTime:departureTime]=0
    Aeq_p = np.diagflat(a.to_numpy())
    Aeq_p = np.concatenate((Aeq_p,np.zeros((length,3*length))),axis=1)
    
    for i in range(len(Aeq_p)):
        if Aeq_p[i,i] ==0 : # delete empty lines
            pass
        else:    
            def feq_j(j):
                def feq(x):
                    return np.dot(Aeq_p[j],x)
                return feq
            eq_constraints.append({'type':'eq', 'fun': feq_j(i)})

        # p+q+lambda+r = 0
    Aeq_sum = np.concatenate((np.eye(length),np.eye(length),np.eye(length),np.eye(length)),axis=1)
    for i in range(len(Aeq_sum)):
        def feq_j(j):
            def feq(x):
                return np.dot(Aeq_sum[j],x)
            return feq
        eq_constraints.append({'type':'eq', 'fun': feq_j(i)})
    
    
    # inequality constraint:
    ineq_constraints = []
    
        # SOCmin <= SOC(t) <= SOCmax
    A1Low = (etaEV/EmaxEV)*(timestep/60)*np.tril(np.ones((length,length)))
    A1Low = np.concatenate((-A1Low,np.zeros((length,3*length))),axis=1)
    A1Up = -A1Low
    B1Low = (SOCmin-SOCinit) * np.ones((length,1))
    B1Up = -(SOCmax-SOCinit) * np.ones((length,1))
    
        # SOCfinal >= SOCrequired
    A2 = -(etaEV/EmaxEV)*(timestep/60)*np.ones((1,length))
    A2 = np.concatenate((A2,np.zeros((1,3*length))),axis=1)
    B2 = (SOCrequired - SOCinit) * np.ones((1,1))

    A = np.concatenate((A1Low,A1Up,A2))
    B = np.concatenate((B1Low, B1Up,B2))
    
    for i in range(len(A)):
        def fineq_j(j):
            def fineq(x):
                return np.dot(A[j],x)-B[j]
            return fineq
        ineq_constraints.append({'type':'ineq', 'fun': fineq_j(i)})

    # initial vector
    x0 = yK 

    
#    res = opt.minimize(utilityEV, x0=x0, method='SLSQP', bounds = bnds,
#                       constraints= ineq_constraints + eq_constraints,
#                       options = {'disp' : True, 'maxiter': 150, 'disp':True}, tol = 0.5)
    
    res = opt.minimize(utilityEV, x0=x0, method='SLSQP', bounds = bnds,
                       constraints= ineq_constraints + eq_constraints,
                       options = {'disp' : True, 'maxiter': 500, 'disp':True}, tol = 0.005)
        
    

    
    yEV = res.x
    return yEV


###############################################################################

def optimStor(beta,parametersSt,startTime,endTime,timestep,price_imp,price_PV,price_B,rho,L,yK, yPVforecast):  
    
    """ optimization of the storage utility function with penalty parameter rho
    
    
        U_storage(x) = - beta1*(pi_local*q+pi_imp*lambda)
                       - beta2*(p_PV_forecast-q)**2
        function to minimize: -sum(U_storage(x)) + (rho/2)*||x-yK+L||**2

        with variables:
            p: net production (<0 if consumption)
            q: energy imported from the energy community (>0 if imported, <0 if exported)
            lambda: energy imported from external grid (>0)
            r: energy imported from the battery (>0 if imported, <0 if exported)
        
        
        Inputs: 
            - beta = [beta1,beta2] (costs, local exchanges)
            - parametersSt = [PmaxIn,PmaxOut,Emax,eta,SOCmin,SOCmax,SOCinit]
            - startTime, endTime: pd.to_datetime()
            - timestep: minutes
            - price_imp, price_PV,  yPVforecast: pd.Series
            - yK, L: np.array
            - rho: floats
        Outputs:
            - yStor: np.array - variables (p,q,lambda) that maximize utitility    
    """    

    PmaxIn = parametersSt[0]
    PmaxOut = parametersSt[1]
    Emax = parametersSt[2]
    eta = parametersSt[3]
    SOCmin = parametersSt[4]
    SOCmax = parametersSt[5]
    SOCinit = parametersSt[6]

    timeVector = pd.date_range(startTime, endTime,freq = str(timestep)+'T')
    length = len(timeVector) # vector length
    
    # selection on the correct time interval
    price_imp = price_imp[startTime:endTime] 
    price_PV = price_PV[startTime:endTime]
    price_B = price_B[startTime:endTime]
    
    yPVforecast =  yPVforecast[startTime:endTime]
    yPVforecast = np.array(yPVforecast)
    yPVforecast = yPVforecast+0.000001  # avoid 0 as it creates problem that I don't understand
        
    # intial step
    if len(yK) !=0: 
        pass
    else:
        #yK=np.zeros(4*length) # initial step
        #yK = 0.0001*np.ones(4*length)
        yK = np.ones(4*length)
        
    # utility function to optimize
    def utilityStor(x):
        """
        U_storage(x) = - beta1*(pi_local*q+pi_imp*lambda)
                       - beta2*(p_PV_forecast-q)**2
        
        function to minimize: -sum(U_storage(x)) + (rho/2)*||x-yK+L||**2
        
        p = x[0:length]
        q = x[length:2*length]
        lambda = x[2*length:3*length]
        r = x[3*length:4*length]
        """
        utility = -beta[0]*(price_PV*x[length:2*length]+price_imp*x[2*length:3*length]+price_B*x[3*length:4*length])/((np.max(np.maximum.reduce([price_imp, price_PV, price_B])))*np.maximum(PmaxIn,PmaxOut)) \
                  -beta[1]*((yPVforecast-x[length:2*length])**2)/np.amax(yPVforecast-np.maximum(PmaxIn,PmaxOut))**2 
                  
        sumUtility = - np.sum(utility)+ (rho/2)*np.linalg.norm(x-yK+L)**2
        return(sumUtility)
        
    
        
    # bounds: 
        # -PmaxIn <= p(t) <= PmaxOut
    LB_p = -PmaxIn*np.ones(length)
    UB_p =  PmaxOut*np.ones(length)
    
        # 0 <= q <=  yPVforecast
    LB_q = np.zeros(length)
    UB_q =  yPVforecast
    
        # 0 <= lambda <= PmaxIn_b
    LB_lambda = np.zeros(length)
    UB_lambda = PmaxIn*np.ones(length)
    
        # -PmaxOut_b <= r <= 0
    LB_r = -PmaxOut*np.ones(length)
    UB_r = np.zeros(length)  
    
    LB = np.concatenate((LB_p,LB_q,LB_lambda,LB_r))
    UB = np.concatenate((UB_p,UB_q,UB_lambda,UB_r))

    bnds = opt.Bounds(LB,UB)    
    
    
    # equality constraint: 
    eq_constraints = []    
    
        # SOCfinal = SOCinit <=> sum(p)=0
    Aeq_p = np.ones((1,length))
    Aeq_p = np.concatenate((Aeq_p,np.zeros((1,3*length))),axis=1)
    eq_constraints.append({'type':'eq', 'fun' : lambda x: np.dot(Aeq_p,x)})
    
    #     # p+q+lambda+r= 0
    Aeq_sum = np.concatenate((np.eye(length),np.eye(length),np.eye(length),np.eye(length)),axis=1)
    for i in range(len(Aeq_sum)):
        def feq_j(j):
            def feq(x):
                return np.dot(Aeq_sum[j],x)
            return feq
        eq_constraints.append({'type':'eq', 'fun': feq_j(i)})
		
	# test avec r**2 + (q+lambda)**2 = p**2
    # for i in range(length):
    #     def feq_nonlin_qlambdar_j(j):
    #         def feq_nonlin_qlambdar(x):
    #             return (-x[j]**2+(x[length+j]+x[2*length+j])**2+x[3*length+j]**2)
    #         return feq_nonlin_qlambdar
    #     eq_constraints.append({'type':'eq', 'fun': feq_nonlin_qlambdar_j(i)})              
    
    
    
    # # inequality constraint: 
    ineq_constraints = []
    
        # SOCmin <= SOC(t) <= SOCmax
    A1Low = -(eta/Emax)*(timestep/60)*np.tril(np.ones((length,length)))
    A1Low = np.concatenate((A1Low,np.zeros((length,3*length))),axis=1)
    A1Up = -A1Low
    B1Low = (SOCmin-SOCinit) * np.ones((length,1))
    B1Up = -(SOCmax-SOCinit) * np.ones((length,1))
    
    A = np.concatenate((A1Low,A1Up))
    B = np.concatenate((B1Low, B1Up))
    
    for i in range(len(A)):
        def fineq_j(j):
            def fineq(x):
                return np.dot(A[j],x)-B[j]
            return fineq
        ineq_constraints.append({'type':'ineq', 'fun': fineq_j(i)})
        
        
        
        # test avec r**2 + (q+lambda)**2 <= p**2
    for i in range(length):
        def feq_nonlin_qlambdar_j(j):
            def feq_nonlin_qlambdar(x):
                return (x[j]**2-(x[length+j]+x[2*length+j])**2-x[3*length+j]**2)
            return feq_nonlin_qlambdar
        ineq_constraints.append({'type':'ineq', 'fun': feq_nonlin_qlambdar_j(i)})                 
        
        
        
        
    #       #  p^2 - (q+lambda)2 > 0 (p^2>(q+lambda)^2)
    # for i in range(length):
    #     def fineq_nonlin_pqlambda_j(j):
    #         def fineq_nonlin_pqlambda(x):
    #             return (x[j]**2-(x[length+j]+x[2*length+j])**2)
    #         return fineq_nonlin_pqlambda
    #     ineq_constraints.append({'type':'ineq', 'fun': fineq_nonlin_pqlambda_j(i)})
            
    #     # p^2 - (r)2 > 0> 0 (p^2>r^2)
    # for i in range(length):
    #     def fineq_nonlin_pr_j(j):
    #         def fineq_nonlin_pr(x):
    #             return (x[j]*2-x[3*length+j]*2)
    #         return fineq_nonlin_pr
    #     ineq_constraints.append({'type':'ineq', 'fun': fineq_nonlin_pr_j(i)})   

    # initial vector
    x0 = np.concatenate((yK[0:length],yPVforecast,yK[2*length:4*length]))

    
    res = opt.minimize(utilityStor, x0=x0, method='SLSQP', bounds = bnds,
                        constraints= ineq_constraints + eq_constraints,
                        options = {'disp' : True, 'maxiter': 500}, tol = 1)

    #res=opt.differential_evolution(utilityStor, bounds = bnds, constraints= ineq_constraints + eq_constraints)
    yStor = res.x
    return yStor   


###############################################################################
    
def optimPV(gamma,startTime,endTime,timestep,p_forecast,price_imp,price_PV,rho,L,yK):
    
    """ optimization of the PV utility function with penalty parameter rho
        U_PV(x) = - gamma_1*(pi_PV*q)
                  - gamma_2*(porecast-q)**2
                  - gamma_3*(p^forecast-p)**2
    
        function to minimize: -sum(U_PV(x)) + (rho/2)*||x-yK+L||**2
       
        with variables:
            p: net production (<0 if consumption)
            q: energy imported from the energy community (>0 if imported, <0 if exported)
            lambda: energy imported from external grid (>0)
            r: energy imported from the battery (>0 if imported, <0 if exported)
        

       note:
       for PV generators, we automatically have lambda = 0 and r = 0 (no energy import)
       thus p = -q

        Inputs: 
            - gamma = [gamma1,gamma2,gamma3] (costs, local exchanges,comfort)
            - startTime, endTime: pd.to_datetime()
            - timestep: minutes
            - p_forecast, price_imp, price_PV: pd.Series
            - Yk, L: np.array
            - rho: floats
        Outputs:
            - yPV: np.array - variables (p,q,lambda) that maximize utitility              
    """    
    
    timeVector = pd.date_range(startTime, endTime,freq = str(timestep)+'T')
    length = len(timeVector) # vector length
    
    # selection on the correct time interval
    p_forecast = p_forecast[startTime:endTime] 
    p_forecast = np.array(p_forecast)

    price_imp = price_imp[startTime:endTime]
    price_PV = price_PV[startTime:endTime]
    
    # intial step
    if len(yK) !=0: 
        pass
    else:
        #yK=np.zeros(4*length) # initial step
        yK = np.concatenate((p_forecast, -p_forecast,np.zeros(2*length)))
        
    p_forecast = p_forecast+0.000001 # avoid 0 as it creates problem that I don't understand

    # utility function to optimize
    def utilityPV(x):
        """
        U_PV(x) = - gamma_1*(pi_PV*q)
                  - gamma_2*(porecast-q)**2
                  - gamma_3*(p^forecast-p)**2
    
        function to minimize: -sum(U_PV(x)) + (rho/2)*||x-yK+L||**2
        
        p = x[0:length]
        q = x[length:2*length]
        lambda = x[2*length:3*length]
        r = x[3*length:4*length]
        """
        #utility = - gamma[0]*(price_PV*x[length:2*length]) \
        utility = - gamma[0]*(price_PV*x[length:2*length])/(np.maximum(price_PV,price_imp)*np.max(p_forecast))\
                  - gamma[1]*((p_forecast-x[length:2*length])**2)/(np.max(p_forecast)**2) \
                  - gamma[2]*((p_forecast-x[0:length])**2)/(np.max(p_forecast))**2
              
        sumUtility = - np.sum(utility)+ (rho/2)*np.linalg.norm(x-yK+L)**2
        return(sumUtility)
    
    # bounds: 
        # 0< p(t) < p_forecast
    LB_p = np.zeros(length)
    UB_p = p_forecast
    
        # q<0 (no import from community)
    LB_q = -np.inf*np.ones(length)    
    UB_q = np.zeros(length)
    
        # lambda = 0, r = 0
    LB_lambda_r = np.zeros(2*length)
    UB_lambda_r = 0.00001*np.ones(2*length)
    
    LB = np.concatenate((LB_p,LB_q,LB_lambda_r))
    UB = np.concatenate((UB_p,UB_q, UB_lambda_r))
    bnds = opt.Bounds(LB,UB)    
    
    # equality constraints:
    eq_constraints = [] 
    
        # lambda = 0

        # p+q= 0
    Aeq_sum = np.concatenate((np.eye(length),np.eye(length),np.zeros((length,length)),np.zeros((length,length))),axis=1)
    for i in range(len(Aeq_sum)):
        def feq_j(j):
            def feq(x):
                return np.dot(Aeq_sum[j],x)
            return feq
        eq_constraints.append({'type':'eq', 'fun': feq_j(i)})
           
    
    
    x0 = np.concatenate((p_forecast,-p_forecast,0.0001*np.zeros(2*length))) # intial vector
    
    res = opt.minimize(utilityPV, x0=x0, method='SLSQP', bounds = bnds,
                       constraints = eq_constraints, options = {'disp' : True, 'maxiter':150}, tol=0.5)

    yPV = res.x
    return yPV


###############################################################################

def optimConsumer(delta,flex,p_forecast,startTime,endTime,timestep,price_imp,price_PV,price_B,rho,L,yK, yPVforecast):
    
    """ optimization of the consumer utility function with penalty parameter rho
        
        U_Load(x)  = - delta_1*(pi_PV*q+pi_imp*lambda+pi_b*r)
                     - delta_2*(p_forecast_PV_total-q)**2
                     + delta_3*(p^forecast-p)**2

        function to minimize: -sum(U_Load(x)) + (rho/2)*||x-yK+L||**2
        
        with variables:
            p: net production (<0 if consumption)
            q: energy imported from the energy community (>0 if imported, <0 if exported)
            lambda: energy imported from external grid (>0)
            r: energy imported from the battery (>0 if imported, <0 if exported)
        

        Inputs: 
            - delta = [delta1,delta2,delta3] (costs, local exchanges,comfort)
            - startTime, endTime: pd.to_datetime()
            - timestep: minutes
            - p_forecast, price_imp, price_PV: pd.Series
            - Yk, L: np.array
            - rho: floats
            -  yPVforecast: np.array, total PV forecast
        Outputs:
            - yLoad: np.array - variables (p,q,lambda) that maximize utitility              
    """    
    
    timeVector = pd.date_range(startTime, endTime,freq = str(timestep)+'T')
    length = len(timeVector) # vector length
    
    # selection on the correct time interval
    p_forecast = -np.abs(p_forecast[startTime:endTime]) 
    p_forecast = np.array(p_forecast)
    price_imp = price_imp[startTime:endTime]
    price_PV = price_PV[startTime:endTime]
    price_B = price_B[startTime:endTime]
    
    yPVforecast = np.abs(yPVforecast[startTime:endTime])
    yPVforecast = np.array(yPVforecast)
    yPVforecast = yPVforecast+0.000001 # avoid 0 as it creates problem that I don't understand
    
    
    # intial step
    if len(yK) !=0: 
        pass
    else:
        #yK=np.zeros(4*length) # initial step  
        #yK = 0.000001*np.ones(4*length)
        yK = np.concatenate((p_forecast,-p_forecast,np.zeros(2*length)))


	# utility function to optimize
    def utilityLoad(x):
        """
        U_Load(x)  = - delta_1*(pi_PV*q+pi_imp*lambda+pi_b*r)
                     - delta_2*(p_forecast_PV_total-q)**2
                     + delta_3*(p^forecast-p)**2

        function to minimize: -sum(U_Load(x)) + (rho/2)*||x-yK+L||**2

        p = x[0:length]
        q = x[length:2*length]
        lambda = x[2*length:3*length]    
        r = x[3*length:4*length] 
        """
        utility = - 0*(price_PV*x[length:2*length]+price_imp*x[2*length:3*length]+price_B*x[3*length:4*length])/np.abs((np.maximum.reduce([price_imp, price_PV, price_B]))*np.max(p_forecast))\
                  - 0*((yPVforecast-x[length:2*length])**2)/(np.max(yPVforecast-p_forecast)**2) \
                  - 1*((p_forecast-x[0:length])**2)/(np.max(p_forecast))**2
        
        sumUtility = - np.sum(utility)+ (rho/2)*np.linalg.norm(x-yK+L)**2
        return(sumUtility)
    
    
    
    
    # bounds: 
        # (1-flex)*p_forecast >= p(t) >= (1+flex)*p_forecast
    LB_p = (1+flex)*p_forecast
    UB_p = (1-flex)*p_forecast
    
        # q>0 (no export to community)
    LB_q = np.zeros(length)
    UB_q =  yPVforecast
    
        # 0 <= lambda 
    LB_lambda = np.zeros(length)
    UB_lambda = np.inf*np.ones(length)
    
            # 0 <= r
    LB_r = np.zeros(length)
    UB_r = np.inf*np.ones(length)
    
    LB = np.concatenate((LB_p,LB_q,LB_lambda,LB_r))
    UB = np.concatenate((UB_p,UB_q,UB_lambda,UB_r))
    bnds = opt.Bounds(LB,UB)    
    
  	# equality constraints:
    eq_constraints = []     
    
        # total energy constant <=> sum(p) = sum(p_forecast)
    Aeq_p = np.ones((1,length))
    Aeq_p = np.concatenate((Aeq_p,np.zeros((1,3*length))),axis=1)
    eq_constraints.append({'type':'eq', 'fun' : lambda x: np.dot(Aeq_p,x)-np.sum(p_forecast)})
    
        # p+q+lambda +r= 0
    Aeq_sum = np.concatenate((np.eye(length),np.eye(length),np.eye(length),np.eye(length)),axis=1)
    for i in range(len(Aeq_sum)):
        def feq_j(j):
            def feq(x):
                return np.dot(Aeq_sum[j],x)
            return feq
        eq_constraints.append({'type':'eq', 'fun': feq_j(i)})
           
    #x0 = np.concatenate((p_forecast,yK[length:4*length])) # intial vector
    x0 = np.concatenate((p_forecast,5*np.ones(3*length)))
						  
    res = opt.minimize(utilityLoad, x0=x0, method='SLSQP', bounds = bnds,
                       constraints = eq_constraints,
                       options = {'disp' : True, 'maxiter' : 500}, tol = 0.05)

    yLoad = res.x
    return yLoad
