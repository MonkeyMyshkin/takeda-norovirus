%Script to build Vaccine ODE function
%this should only be needed once in order to build the ode function for
%simulation 
clear all
Lmax=81;


%DECLARE ageGroupBreaks
ageGroupBreaks=[8 18 26 37 50 70];
noAgeGroups=length(ageGroupBreaks)+1;

%DECLARE symbolic variables
syms t theta
StratifiedContactMatrix =sym('StratifiedContactMatrix',[length(ageGroupBreaks)+1,length(ageGroupBreaks)+1]);
mu=sym('mu',[1,Lmax]);
phi=sym('phi',[1,Lmax]);
param=sym('param',[1,15]);
omega2=sym('omega2',[1,7]);

%CALL MakeVaccineEquations function to create symbolic equations
[Equations] = MakeEquationsVaccination(mu,param(1),param(2),param(3), param(4), param(5), param(6), param(7), param(8), ...
    param(9),param(10) ,param(11) ,param(12) ,param(13) ,param(14) ,omega2 ,phi , theta, StratifiedContactMatrix ,ageGroupBreaks);

%REPLACE symbolic variable names with one variable, Y
M=sym('M',[Lmax,1]);S=sym('S',[Lmax,1]);E1=sym('E1',[Lmax,1]);E2=sym('E2',[Lmax,1]);I=sym('I',[Lmax,1]);A=sym('A',[Lmax,1]);R=sym('R',[Lmax,1]);
Sv=sym('Sv',[Lmax,1]);E1v=sym('E1v',[Lmax,1]);E2v=sym('E2v',[Lmax,1]);Iv=sym('Iv',[Lmax,1]);Av=sym('Av',[Lmax,1]);Rv=sym('Rv',[Lmax,1]);
%ode variable
Y=sym('Y',[13*Lmax,1]);


%sub in ode variable
EQNS=subs(Equations,[M.',S.',E1.',E2.',I.',A.',R.',Sv.',E1v.',E2v.',Iv.',Av.',Rv.'],Y.');


%use MATLAB built in function to optimise model equations for use with ODE
%solvers
matlabFunction(EQNS, 'file', 'NoroModelVaccODE','vars', {t, Y,param,omega2,mu ,phi,theta,StratifiedContactMatrix});
disp('Built ODE file')