function [ TIME,SimulationResult ] = SimulateSeasons_Vaccination( Param,ParamVac,omega2,mu ,phi,theta, NewContactMatrix , x0)
%Take input parameters and initial conditions and simulate model for all
%seasons and all age groups

Lmax=length(x0)/7;

%transform so gamma is scaling of duration of symptomatic duration
Param(9)=Param(8)*Param(9);

options=odeset('NonNegative',1:15*Lmax);%,'RelTol',1e-5,'AbsTol',1e-5);

Parameters=[Param, ParamVac];
%measure one season from week 33 to week 33 of following year
%time increments are 7 days
CM=reshape(NewContactMatrix,[],1)/1e6;

%SEASON ONE
tspan=1:1:33*7;
[TIME,X]=ode45(@(t,Y)NoroModelAR_Vacc_ODE(t,Y,Parameters ,omega2,mu,phi,theta(1),CM),tspan,[1e6*x0,zeros(1,9*Lmax)],options);

%SEASON TWO to EIGHT

for seasonIndex= 2 : 8
    
    %simulate for each season
    
    [ TIMEnew , Xnew ]=ode45(@(t,Y)NoroModelAR_Vacc_ODE(t,Y,Parameters,omega2,mu ,phi,theta(seasonIndex),CM),(TIME(end):1:TIME(end)+52*7),X(end,:),options);
    
    X = [ X ; Xnew(2:end,:) ];
    TIME = [ TIME ; TIMEnew(2:end) ];
end

%SEASON NINE
[ TIMEnew , Xnew ]=ode45(@(t,Y)NoroModelAR_Vacc_ODE(t,Y,Parameters ,omega2,mu,phi,theta(9), CM),(TIME(end):1:TIME(end)+21*7),X(end,:),options);
SimulationResult = [ X ; Xnew(2:end,:) ];
TIME = [ TIME ; TIMEnew(2:end) ];


    
maxInd=length(TIME);

 TIME=TIME(2:7:maxInd);
 SimulationResult=SimulationResult(2:7:maxInd,:)/1e6;
%  
 
 
end

