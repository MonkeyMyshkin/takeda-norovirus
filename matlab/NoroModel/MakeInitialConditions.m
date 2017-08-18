function [ InitialConditions] = MakeInitialConditions( Param,theta,mu,ContactMatrix )
%%
%MakeInitialConditions simulates models to arrive at stable age
%distribution, creating x0 for simulateSeasons

%INPUTS
%Param=vector of parameters, minimum requirement is alpha, q, omega1,
%      omega2, nu, delta, epsilon, sigma, psi, gamma. 
%theta=susceptibilities of recoveed individuals per season
%mu= death rates per age
%ContactMatrix= probabilities of contact between age groups

%OUTPUTS
%InitialConditions=x0 for SimulateSeasons
%%
Lmax=length(mu);

load('fixedPoints.mat') %Output from previous simulation as startpoint

options=odeset('NonNegative',1:9*Lmax);%ODE solver option

%contact matrix entered as vector
CM=reshape(ContactMatrix,[],1)/1e6;

%Burn in period of 50 years
tspan=1:1:50*365;

[TIME,X]=ode45(@(t,Y)NoroModelARODE(t,Y,Param,mu,theta(1),CM),tspan,[1e6*fixedPoints,zeros(1,2*Lmax)],options);  %ICSarescaled up to improve accuracy

%time to average of 2 years

[TIMENEW , Xnew ]=ode45(@(t,Y)NoroModelARODE(t,Y,Param,mu ,theta(1),CM),(TIME(end):1:TIME(end)+2*365),X(end,:),options);

%
InitialConditions = mean(Xnew(:,1:7*Lmax));

InitialConditions=InitialConditions/1e6;    %scale output back down
end

