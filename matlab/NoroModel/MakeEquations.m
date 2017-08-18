function [ Equations ] = MakeEquations(param,mu,theta,StratifiedContactMatrix,ageGroupBreaks )
%%
%MakeEquations makes the equations to be build into an ODE script
%INPUTS
%param= parameter vector in following order
%   alpha=loss of maternal antibodies
%   q=tranmissability
%   omega1= seasonal amplitude
%   omega2=seasonal offset
%   nu=scaling of asymptomatic infectiousness
%   delta=loss of immunity
%   epsilon=rate of latency loss
%   sigma=proportion symptomatic
%   psi=rate infected individuals become asymptomatic
%   theta= susceptibility
%   chi= rate innate immunity is lost
%mu=death rate
%StratifiedContactMatrix= Mixing between age groups, split by age group not
%                         year
%ageGroupBreaks=a vector of age group divisions 

%OUTPUTS
%Equations=System of differential equations

%%
%assign parameter names
alpha=param(1); q=param(2); omega1=param(3); omega2=param(4); nu=param(5); 
delta=param(6); epsilon=param(7); sigma=param(8); psi=param(9); 
gamma=param(10); 

C=StratifiedContactMatrix;


%define maximum age
Lmax=length(mu);
ageGroupBreaks=[ageGroupBreaks,Lmax];
NoAgeGroups=length(ageGroupBreaks);

%age rate for yearly ageing in terms of days
age=1/365;  

%DECLARE symbolic variables
M=sym('M',[Lmax,1]);S=sym('S',[Lmax,1]);E1=sym('E1',[Lmax,1]);E2=sym('E2',[Lmax,1]);I=sym('I',[Lmax,1]);A=sym('A',[Lmax,1]);R=sym('R',[Lmax,1]);Im=sym('Im',[Lmax,1]);
syms t


%seasonal forcing term
Z=(1+omega1*cos(2*pi*t/364+omega2));

%calculate force of infection
for index=1:Lmax
    %find which age group index is in
    ageGroupIndex=find(histcounts(index,[-1,ageGroupBreaks]));
    %force of infection depends on transmission parameter q, seasonal
    %forcing term Z, contact matrix C, infectious symptomatic class I and
    %infectious asymptomatic class A
    FOI(index)= q*Z*C(ageGroupIndex,1)*sum(I(1:ageGroupBreaks(1))+nu*A(1:ageGroupBreaks(1)));
    for jndex = 2 : NoAgeGroups
        FOI(index)=FOI(index) ...
                    + q*Z*C(ageGroupIndex,jndex)*sum(I(ageGroupBreaks(jndex-1)+1:ageGroupBreaks(jndex))+nu*A(ageGroupBreaks(jndex-1)+1:ageGroupBreaks(jndex)));
    end
end

%closed population so births equal all deaths
BIRTHS=sum(mu.*(M+S+E1+E2+I+A+R).');

%EQUATIONS
%age 1
mat(1)=BIRTHS - (alpha +mu(1))*M(1) - age*M(1) ;
sus(1)=alpha*M(1) - (FOI(1)+mu(1))* S(1) + delta*R(1) - age*S(1);
exposed1(1)= FOI(1)*S(1) - (epsilon+mu(1))*E1(1) - age*E1(1);
exposed2(1)= epsilon*E1(1) - (epsilon+mu(1))*E2(1) - age*E2(1);
inf(1)= sigma*epsilon*E2(1) - (psi+mu(1))*I(1) - age*I(1);
asymp(1)= (1-sigma)*epsilon*E2(1) + psi*I(1) + FOI(1)*theta*R(1) - (gamma+mu(1))*A(1) - age*A(1);
rec(1)= gamma*A(1) - FOI(1)*theta*R(1) - (delta+mu(1))*R(1) - age*R(1);
%bookkeeping, efflux from infectious class
Cas(1)=psi*I(1);
AllCases(1)=epsilon*E2(1);

%age 2 and above
for index=2:Lmax-1

    mat(index)= -(alpha +mu(index))*M(index) + age*(M(index-1)-M(index)) ;
    
    sus(index)=alpha*M(index) - (FOI(index)+mu(index))* S(index) + delta*R(index) + age*(S(index-1)-S(index));
    
    exposed1(index)= FOI(index)*S(index) - (epsilon+mu(index))*E1(index) + age*(E1(index-1)-E1(index));
    
    exposed2(index)= epsilon*E1(index) - (epsilon+mu(index))*E2(index) + age*(E2(index-1)-E2(index));
    
    inf(index)= sigma*epsilon*E2(index) - (psi+mu(index))*I(index) + age*(I(index-1)-I(index));
    
    asymp(index)= (1-sigma)*epsilon*E2(index) + psi*I(index) + FOI(index)*theta*R(index) - (gamma+mu(index))*A(index) + age*(A(index-1)-A(index));
    
    rec(index)=gamma*A(index) - FOI(index)*theta*R(index) - (delta+mu(index))*R(index) + age*(R(index-1)-R(index));
    
    Cas(index) = psi*I(index);
    
    AllCases(index)=epsilon*E2(index);
end

%last age 
mat(Lmax)= - (alpha +mu(Lmax))*M(Lmax) + age*(M(Lmax-1)) ;
sus(Lmax)=alpha*M(Lmax) - (FOI(Lmax)+mu(Lmax))* S(Lmax) + delta*R(Lmax) + age*(S(Lmax-1));
exposed1(Lmax)= FOI(Lmax)*S(Lmax) - (epsilon+mu(Lmax))*E1(Lmax) + age*(E1(Lmax-1));
exposed2(Lmax)= epsilon*E1(Lmax) - (epsilon+mu(Lmax))*E2(Lmax) + age*(E2(Lmax-1));
inf(Lmax)= sigma*epsilon*E2(Lmax) - (psi+mu(Lmax))*I(Lmax) + age*(I(Lmax-1));
asymp(Lmax)= (1-sigma)*epsilon*E2(Lmax) + psi*I(Lmax) + FOI(Lmax)*theta*R(Lmax) - (gamma+mu(Lmax))*A(Lmax) + age*(A(Lmax-1));
rec(Lmax)=gamma*A(Lmax) - FOI(Lmax)*theta*R(Lmax) - (delta+mu(Lmax))*R(Lmax) + age*(R(Lmax-1));
Cas(Lmax) = psi*I(Lmax);
AllCases(Lmax)=epsilon*E2(Lmax);

%concatenate equations
Equations=[mat.';sus.';exposed1.';exposed2.';inf.';asymp.';rec.';Cas.';AllCases.'];
end

