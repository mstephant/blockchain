%% Battery local optimization for ADMM
% optimization according to the utility function of the battery
% utility function is : 
% U = -beta1*(pi_PV*q+pi_imp*lambda-pi_b*r)
%     -beta2*(p_PV_forecast-q^2)

% with variables : 
%   p: net production (<0 if consumption)
%   q: energy imported from the energy community (>0 if imported, <0 if
%   exported)
%   lambda: energy imported from external grid (>0)
%   r: energy imported from the battery (>0 if imported, <0 if exported)

function [X_B,utility_value] = util_func_B_4(beta1,beta2,parametersB,timestep,startTime,endTime,price_PV,price_b,price_imp,rho,L,X_k,p_forecast_PV_total)

    % 
    % INPUTS :
    %   - beta1, beta2 (real numbers in [0,1]) : coefficients de préférence
    %   (coût, échanges avec la communauté)
    %   - parametersB = [PmaxIn_b,PmaxOut_b,Emax_b,eta_b,SOCmin_b,SOCmax_b,SOCinit_b]
    %                   paramètres de la batterie
    %   - timestep (duration) : pas de temps des vecteurs
    %   - startTime (datetime) : instant du début d'optimisation
    %   - endTime (datetime) : instant de fin d'optimisation
    %   - price_PV (timetable) : prix de l'électricité PV (€/kWh)
    %   - price_b (timetable) : prix de l'électricité de la batterie(€/kWh) 
    %   - price_imp (timetable) : prix de l'électricité importée du
    %                             réseau externe (€/kWh)
    %   - L: value of mean(x^k) - mean(z^k) + u^k (from previous iteration of ADMM)
    %   - rho: stpe of ADMM algorithm 
    %   - X_k: solution of the previous ADMM iteration for player i only
    %   - p_forecast_PV_total (table)
    %
    % OUTPUTS :
    %   - X_B : paramètres (p,q,lambda) qui maximisent l'utilité (kW)
    %   - utility_B : valeur de l'utilité (satisfaction de l'usager)
    
    
    %% Définition des paramètres 
    PmaxIn_b = parametersB(1); % (kW) puissance max en charge de la batterie
    PmaxOut_b = parametersB(2); % (kW) puissance min en décharge de la batterie (doit être positive)
    Emax_b = parametersB(3); % (kWh) capacité de la batterie
    eta_b = parametersB(4); % rendement de la batterie
    SOCmin_b = parametersB(5); % etat de charge minimal
    SOCmax_b = parametersB(6); % etat de charge maximal
    SOCinit_b = parametersB(7); % etat de charge initial
    
    timeVector_b = startTime:timestep:endTime;   
    timeRange_b = timerange(startTime,endTime,'closed');
    price_PV = price_PV(timeRange_b,:);% sélection sur l'intervalle de temps
    price_imp = price_imp(timeRange_b,:);
    price_b = price_b(timeRange_b,:);

    N = length(timeVector_b);
    
    %% Définition des contraintes
    
    % Bornes
    % -PmaxIn_b <= p(t) <= PmaxOut_b
    LB_p = -PmaxIn_b * ones(N,1);  % limite inférieure pour la puissance 
    UB_p = PmaxOut_b * ones(N,1);   % limite supérieure pour la puissance


     % la contrainte ci-dessous sont redondantes avec la première :
 %     % 0 <= q <= p_forecast_PV_total
     LB_q = zeros(N,1);
     UB_q = abs(p_forecast_PV_total);
     
     % 0<= lambda <= PmaxIn_b
     LB_lambda = zeros(N,1);
     UB_lambda = PmaxIn_b*ones(N,1);

     % -PmaxOut_b <= r <= 0
     LB_r = -PmaxOut_b * ones(N,1);
     UB_r = zeros(N,1);

     LB = [LB_p;LB_q;LB_lambda;LB_r];
     UB = [UB_p;UB_q;UB_lambda;UB_r]; 

    % Contraintes d'inégalité
    % SOC borné
    % SOC(t+1) = SOCinit_b + eta_b/E_b *sum(power(1:t))*t
    b = [(SOCmax_b - SOCinit_b)*ones(N,1); (-SOCmin_b + SOCinit_b)*ones(N,1)];
    s = eta_b/Emax_b*hours(timestep)*sparse(tril(ones(N))) ;
    s = [s,zeros(N,3*N)];
    A = [-s;s];
    
    % Contraintes d'égalité :
    % SOC(final) = SOC(initial) ==> Sum(p_i)=0
    Aeq_p = ones(1,N);
    Aeq_p = [Aeq_p,zeros(1,3*N)];
    beq_p = zeros(1,1);
    
    % contrainte: p+q+lambda+r = 0
    Aeq_sum = [eye(N),eye(N),eye(N),eye(N)];
    beq_sum = zeros(N,1);
    
    % concaténation
    Aeq = [Aeq_p;Aeq_sum];
    beq = [beq_p;beq_sum];
    
    % contrainte non linéaire (pour éviter spéculation) : 
%     % - p^2+q^2 <= 0
%     % - p^2+lambda^2 <= 0 
%     c = @(x)[-x(1:N).^2+x(N+1:2*N).^2;
%         -x(1:N).^2+x(2*N+1:3*N).^2];

    % q+lambda=0 or q=0 <=> r*(lambda+q)=0
    c = @(x)[];
    ceq = @(x) (x(N+1:2*N)+x(2*N+1:3*N)).*x(3*N+1:4*N);
    %ceq = @(x)[];
    nonlinfcn = @(x)deal(c(x),ceq(x));
     % nonlinfcn = @constraint_storage;
  %% Optimization
  % Fonction objectif:
    %U = -beta1*(pi_PV*q+pi_imp*lambda)-beta2*(p_forecast_PV_total-q)^2    

    % vecteur intial pour l'optimisation
    x0 = X_k;
    
    % Optimize with FMINCON
    options = optimset('MaxFunEvals',Inf,'MaxIter',250,...
    'Algorithm','sqp','Display','iter');

    %options = optimoptions(@fmincon,'Algorithm','sqp-legacy','MaxIterations',500);

    % Rappel : 
    % p = x(1:N);
    % q = x(N+1:2*N);
    % lambda = x(2*N+1:3*N);
    % r = x(3*N+1:4*N);
    
    [X_B,utility_value] = fmincon(@(x) (-sum(-beta1*(price_PV{:,1}.*x(N+1:2*N)...
        +price_imp{:,1}.*x(2*N+1:3*N)+price_b{:,1}.*x(3*N+1:4*N))...
        /(max([price_PV{:,1};price_imp{:,1};price_b{:,1}])*max(PmaxOut_b,PmaxIn_b))...
        -beta2*(p_forecast_PV_total-x(N+1:2*N)).^2/(max(p_forecast_PV_total-max(PmaxOut_b,PmaxIn_b))).^2 ) ...
        + (rho/2)*norm(x-X_k+L)^2),x0,A,b,Aeq,beq,LB,UB,nonlinfcn,options);
    
end 
