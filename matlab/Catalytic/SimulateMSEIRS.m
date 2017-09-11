function [ seropositive,x0] = SimulateMSEIRS( Param, Lmax )
%function to simulate MSEIRS for serological data
%%
%INPUTS
%param          vector of q, delta and nu
%Lmax           maximum age

%OUTPUTS
%seropositive   proportion of each agegroup that are seropositive
%x0             mean final conditions

%%
q=Param(1); delta=Param(2); nu=Param(3);


tspan=0:1:1*365;
options = odeset('NonNegative',(1:5*Lmax));
%initial conditions
load('x0_Son_homog.mat')
x0=x0_Son_homog;
% load('x0_It_homog.mat')
% x0=x0_It_homog;

max_it=10;  %number of years to simulate from initial conditions

for iteration=1:max_it
    
    %solve for one year
    [TIME,x]=ode45(@(t,Y)ODEfile(t,Y,q,delta,nu),tspan,x0,options);
    
    %update all age groups with annual cohort ageing
    BIRTHS=sum(x(end,5*Lmax+1:end));
    
    x0(1:Lmax)= [BIRTHS, x(end,1:Lmax-2),x(end,Lmax-1)+x(end,Lmax)];
    x0(Lmax+1:2*Lmax)= [0, x(end,Lmax+1:2*Lmax-2),x(end,2*Lmax-1)+x(end,2*Lmax)];
    x0(2*Lmax+1:3*Lmax)= [0, x(end,2*Lmax+1:3*Lmax-2),x(end,3*Lmax-1)+x(end,3*Lmax)];
    x0(3*Lmax+1:4*Lmax)= [0, x(end,3*Lmax+1:4*Lmax-2),x(end,4*Lmax-1)+x(end,4*Lmax)];
    x0(4*Lmax+1:5*Lmax)= [0, x(end,4*Lmax+1:5*Lmax-2),x(end,5*Lmax-1)+x(end,5*Lmax)];
    
    %update time
    tspan=(TIME(end):1:TIME(end)+365);
    
end

for iteration=1 %solve for an addition year to calculate average from
    
    %solve
    [TIME1,x1]=ode45(@(t,Y)ODEfile(t,Y,q,delta,nu),tspan,x0,options);
    
    %update all age groups with annual cohort ageing
    BIRTHS=sum(x1(end,5*Lmax+1:end));
    
    x0(1:Lmax)= [BIRTHS, x1(end,1:Lmax-2),x1(end,Lmax-1)+x1(end,Lmax)];
    x0(Lmax+1:2*Lmax)= [0, x1(end,Lmax+1:2*Lmax-2),x1(end,2*Lmax-1)+x1(end,2*Lmax)];
    x0(2*Lmax+1:3*Lmax)= [0, x1(end,2*Lmax+1:3*Lmax-2),x1(end,3*Lmax-1)+x1(end,3*Lmax)];
    x0(3*Lmax+1:4*Lmax)= [0, x1(end,3*Lmax+1:4*Lmax-2),x1(end,4*Lmax-1)+x1(end,4*Lmax)];
    x0(4*Lmax+1:5*Lmax)= [0, x1(end,4*Lmax+1:5*Lmax-2),x1(end,5*Lmax-1)+x1(end,5*Lmax)];
    
    
    tspan=(TIME1(end):1:TIME1(end)+365);
    
end
%Calculate the proportion of each age group that is seropositive
parfor i=1:Lmax
    seropositive(i)=(mean(x1(:,i))+mean(x1(:,2*Lmax+i))+mean(x1(:,3*Lmax+i))+mean(x1(:,4*Lmax+i))+mean(x1(:,5*Lmax+i)))/(mean(x1(:,i))+mean(x1(:,1*Lmax+i))+mean(x1(:,2*Lmax+i))+mean(x1(:,3*Lmax+i))+mean(x1(:,4*Lmax+i))+mean(x1(:,5*Lmax+i)));

end

end

