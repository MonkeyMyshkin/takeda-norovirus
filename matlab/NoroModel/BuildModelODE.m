%Script to build ODE function
%This is used once to create an optimised system of differential equations
%for the model

%start with a clear workspace and DECLARE max age
clear all
Lmax=81;

%DECLARE ageGroupBreaks
ageGroupBreaks=[8 18 26 37 50 70];
noAgeGroups=length(ageGroupBreaks)+1;

%DECLARE symbolic variables
syms t alpha  q omega1 omega2  nu  delta  epsilon  sigma  psi gamma theta
StratifiedContactMatrix =sym('StratifiedContactMatrix',[length(ageGroupBreaks)+1,length(ageGroupBreaks)+1]);
mu=sym('mu',[1,Lmax]);
param=sym('param',[1,10]);

%CALL MakeEquations function to create symbolic equations
[Equations] = MakeEquations(param ,mu, theta, StratifiedContactMatrix ,ageGroupBreaks);

%REPLACE symbolic variable names with one variable, Y
M=sym('M',[Lmax,1]);S=sym('S',[Lmax,1]);E1=sym('E1',[Lmax,1]);E2=sym('E2',[Lmax,1]);I=sym('I',[Lmax,1]);A=sym('A',[Lmax,1]);R=sym('R',[Lmax,1]);
%ode variable
Y=sym('Y',[7*Lmax,1]);
%substitute
EQNS=subs(Equations,[M.',S.',E1.',E2.',I.',A.',R.'],Y.');


%use MATLAB built in function to optimise model equations for use with ODE
%solvers
matlabFunction(EQNS, 'file', 'NoroModelARODE','vars', {t, Y,param,mu ,theta,StratifiedContactMatrix});
