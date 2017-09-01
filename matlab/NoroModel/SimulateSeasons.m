function [ TIME,SimulationResult ] = SimulateSeasons( Param,omega2,mu ,theta, ContactMatrix , x0)
%%
%SimulatSeasons: Takes input parameters and initial conditions and solves the system of ODEs for all
%seasons and all age groups

%INPUTS:
%Param=vector of parameter values. Doubles.
%mu=vector of death rates for each age group. Doubles.
%theta=vector of susceptibility of recovered individuals for each season.
%       Doubles.
%ContactMatrix= noAgeGroups x noAgeGroups matrix of contact probabilities
%               between age groups
%x0=vector of initial conditions for each disease state
%omega2= vector of seasonal offsets for each age group

%OUTPUS:
%TIME=vector of weekly time points
%SimulationResult=solution of ODEs for each week in TIME and each age
%%
Lmax=81;

%options for solver, nonnegative and high degree of accuracy
options=odeset('NonNegative',1:9*Lmax,'RelTol',1e-10,'AbsTol',1e-10);

%measure one season from week 33 to week 33 of following year
%time increments are 7 days
CM=reshape(ContactMatrix,[],1); %matrix enters optimised ode function as vector

%SEASON ONE
tspan=1:1:33*7;
[TIME,X]=ode45(@(t,Y)NoroModelARODE(t,Y,Param,omega2,mu,theta(1),CM),tspan,[x0,zeros(1,2*Lmax)],options);

%SEASON TWO to EIGHT
for seasonIndex= 2 : 8
    
    %simulate for each season
    
    [ TIMEnew , Xnew ]=ode45(@(t,Y)NoroModelARODE(t,Y,Param,omega2,mu,theta(seasonIndex),CM),(TIME(end):1:TIME(end)+52*7),X(end,:),options);
    
    X = [ X ; Xnew(2:end,:) ];
    TIME = [ TIME ; TIMEnew(2:end) ];
end

%SEASON NINE
[ TIMEnew , Xnew ]=ode45(@(t,Y)NoroModelARODE(t,Y,Param,omega2,mu,theta(9),CM),(TIME(end):1:TIME(end)+21*7),X(end,:),options);

%concatenate results
SimulationResult = [ X ; Xnew(2:end,:) ];
TIME = [ TIME ; TIMEnew(2:end) ];

maxInd=length(TIME);

%pick weekly values
TIME=TIME(2:7:maxInd);
SimulationResult=SimulationResult(2:7:maxInd,:);

end

