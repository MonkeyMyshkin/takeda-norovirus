function [ Equations ] = MakeEquationsVaccination(mu, alpha,q,omega1,  nu, delta, epsilon, sigma, psi, gamma, rho, kappa, tau, zeta,sigmaV,omega2,phi,theta, StratifiedContactMatrix ,ageGroupBreaks )
%%
%MakeVaccineEquations makes the vaccine equations to be build into an ODE script


%INPUTS
%Parameter vector- to be ammended
%Contact matrix- split by age group- not annual
C=StratifiedContactMatrix;

%assign parameters from parameter vector!
%alpha=loss of maternal antibodies
%mu=death rate
%q=tranmissability
%nu=scaling of asymptomatic infectiousness
%delta=loss of immunity
%epsilon=rate of latency loss
%sigma=proportion symptomatic
%psi=rate infected individuals become asymptomatic
%omega1, omega2=seasonal forcing terms 
%theta= susceptibility
%phi- vaccination rate including coverage and efficacy
%rho- rate vaccination wanes Sv->S
%kappa- switch between immunizing vaccination and nonimmunizing vaccination
%tau- scaling factor for force of infection experienced by vaccinated
    %individuals
%zeta- scaling of infectiousness of vaccinated individuals
%sigmaV- proportion of vaccination individuals who become symptomatic
%omega2=vector of seasonal offsets

%OUTPUTS
%Equations= system of differential equations for vaccination model

%%
Lmax=length(mu);    %number of age years

ageGroupBreaks=[ageGroupBreaks,Lmax];
NoAgeGroups=length(ageGroupBreaks);

age=1/365;  %age rate for yearly ageing in terms of days


%symbolic variables
M=sym('M',[Lmax,1]);S=sym('S',[Lmax,1]);E1=sym('E1',[Lmax,1]);
E2=sym('E2',[Lmax,1]);I=sym('I',[Lmax,1]);A=sym('A',[Lmax,1]);
R=sym('R',[Lmax,1]);
Sv=sym('Sv',[Lmax,1]);E1v=sym('E1v',[Lmax,1]);E2v=sym('E2v',[Lmax,1]);
Iv=sym('Iv',[Lmax,1]);Av=sym('Av',[Lmax,1]);Rv=sym('Rv',[Lmax,1]);

syms t



Z=(1+omega1.*cos(2*pi*t/364+omega2));

for findex=1:Lmax
    %find which age group index is in
    ageGroupIndex=find(histcounts(findex,[-1,ageGroupBreaks]));
    
    FOI(findex)= q*Z(ageGroupIndex)*C(ageGroupIndex,1)*sum(I(1:ageGroupBreaks(1))+nu*A(1:ageGroupBreaks(1))...
                                    + zeta*(Iv(1:ageGroupBreaks(1))+nu*Av(1:ageGroupBreaks(1))));
    for jndex = 2 : NoAgeGroups
        FOI(findex)=FOI(findex) ...
                    + q*Z(ageGroupIndex)*C(ageGroupIndex,jndex)*sum(I(ageGroupBreaks(jndex-1)+1:ageGroupBreaks(jndex))...
                                                        +nu*A(ageGroupBreaks(jndex-1)+1:ageGroupBreaks(jndex))...
                                                   + zeta*(Iv(ageGroupBreaks(jndex-1)+1:ageGroupBreaks(jndex))...
                                                            +nu*Av(ageGroupBreaks(jndex-1)+1:ageGroupBreaks(jndex))));
    end
end

BIRTHS=sum(mu.*(M+S+E1+E2+I+A+R + Sv+E1v+E2v+Iv+Av+Rv).');

%age 1
mat(1)=BIRTHS - (alpha +mu(1))*M(1) - age*M(1) ;
sus(1)=alpha*M(1) - (FOI(1)+mu(1))* S(1) + delta*R(1) - age*S(1) - phi(1)*S(1) + rho*Sv(1);
exposed1(1)= FOI(1)*S(1) - (epsilon+mu(1))*E1(1) - age*E1(1) - phi(1)*E1(1);
exposed2(1)= epsilon*E1(1) - (epsilon+mu(1))*E2(1) - age*E2(1) - phi(1)*E2(1);
inf(1)= sigma*epsilon*E2(1) - (psi+mu(1))*I(1) - age*I(1);
asymp(1)= (1-sigma)*epsilon*E2(1) + psi*I(1) + FOI(1)*theta*R(1) - (gamma+mu(1))*A(1) - age*A(1) - phi(1)*A(1);
rec(1)= gamma*A(1) - FOI(1)*theta*R(1) - (delta+mu(1))*R(1) - age*R(1) - phi(1)*R(1)+ rho*Rv(1);


susv(1)= kappa*phi(1)*S(1)- (tau*FOI(1)+mu(1))* Sv(1) + delta*Rv(1) - age*Sv(1) - rho*Sv(1);
exposed1v(1)=kappa*phi(1)*E1(1) + tau*FOI(1)*Sv(1) - (epsilon+mu(1))*E1v(1) - age*E1v(1);
exposed2v(1)=kappa*phi(1)*E2(1) + epsilon*E1v(1) - (epsilon+mu(1))*E2v(1) - age*E2v(1);
infv(1)= sigmaV*epsilon*E2v(1) - (psi+mu(1))*Iv(1) - age*Iv(1);
asympv(1)= kappa*phi(1)*A(1) +(1-sigmaV)*epsilon*E2v(1) + psi*Iv(1) + tau*FOI(1)*theta*Rv(1) - (gamma+mu(1))*Av(1) - age*Av(1);
recv(1)= kappa*phi(1)*R(1) +(1-kappa)*phi(1)*(S(1)+E1(1)+E2(1)+A(1)+R(1))+gamma*Av(1) - tau*FOI(1)*theta*Rv(1) - (delta+mu(1))*Rv(1) - age*Rv(1)- rho*Rv(1);
%bookkeeping, efflux from infectious class
Cas(1)=psi*I(1)+psi*Iv(1);
AllCases(1)=epsilon*(E2(1)+E2v(1));

%age 2 and above
for index=2:Lmax-1

    
    mat(index)= -(alpha +mu(index))*M(index) + age*(M(index-1)-M(index)) ;
    
    sus(index)=alpha*M(index) - (FOI(index)+mu(index))* S(index) + delta*R(index) + age*(S(index-1)-S(index)) - phi(index)*S(index)+ rho*Sv(index);
    
    exposed1(index)= FOI(index)*S(index) - (epsilon+mu(index))*E1(index) + age*(E1(index-1)-E1(index))- phi(index)*E1(index);
    
    exposed2(index)= epsilon*E1(index) - (epsilon+mu(index))*E2(index) + age*(E2(index-1)-E2(index))- phi(index)*E2(index);
    
    inf(index)= sigma*epsilon*E2(index) - (psi+mu(index))*I(index) + age*(I(index-1)-I(index));
    
    asymp(index)= (1-sigma)*epsilon*E2(index) + psi*I(index) + FOI(index)*theta*R(index) - (gamma+mu(index))*A(index) + age*(A(index-1)-A(index))...
                    - phi(index)*A(index);
    
    rec(index)=gamma*A(index) - FOI(index)*theta*R(index) - (delta+mu(index))*R(index) + age*(R(index-1)-R(index))- phi(index)*R(index)+ rho*Rv(index);
    
    %vaccinated
    susv(index)=kappa*phi(index)*S(index) - (tau*FOI(index)+mu(index))* Sv(index) + delta*Rv(index) + age*(Sv(index-1)-Sv(index)) - rho*Sv(index);
    
    exposed1v(index)= kappa*phi(index)*E1(index) + tau*FOI(index)*Sv(index) - (epsilon+mu(index))*E1v(index) + age*(E1v(index-1)-E1v(index));
    
    exposed2v(index)= kappa*phi(index)*E2(index) + epsilon*E1v(index) - (epsilon+mu(index))*E2v(index) + age*(E2v(index-1)-E2v(index));
    
    infv(index)= sigmaV*epsilon*E2v(index) - (psi+mu(index))*Iv(index) + age*(Iv(index-1)-Iv(index));
    
    asympv(index)= kappa*phi(index)*A(index) + (1-sigmaV)*epsilon*E2v(index) + psi*Iv(index) + tau*FOI(index)*theta*Rv(index) ...
                - (gamma+mu(index))*Av(index) + age*(Av(index-1)-Av(index));
    
    recv(index)=kappa*phi(index)*R(index) +(1-kappa)*phi(index)*(S(index)+E1(index)+E2(index)+A(index)+R(index))...
                + gamma*Av(index) - tau*FOI(index)*theta*Rv(index) - (delta+mu(index))*Rv(index) + age*(Rv(index-1)-Rv(index)) - rho*Rv(index);
    
    %bookkeeping
    Cas(index) = psi*I(index)+psi*Iv(index);
    
    AllCases(index)=epsilon*(E2(index)+E2v(index));
end

%last age 
mat(Lmax)= - (alpha +mu(Lmax))*M(Lmax) + age*(M(Lmax-1)) ;
sus(Lmax)=alpha*M(Lmax) - (FOI(Lmax)+mu(Lmax))* S(Lmax) + delta*R(Lmax) + age*(S(Lmax-1)) - phi(Lmax)*S(Lmax) +rho*Sv(Lmax);
exposed1(Lmax)= FOI(Lmax)*S(Lmax) - (epsilon+mu(Lmax))*E1(Lmax) + age*(E1(Lmax-1)) -phi(Lmax)*E1(Lmax);
exposed2(Lmax)= epsilon*E1(Lmax) - (epsilon+mu(Lmax))*E2(Lmax) + age*(E2(Lmax-1))-phi(Lmax)*E2(Lmax);
inf(Lmax)= sigma*epsilon*E2(Lmax) - (psi+mu(Lmax))*I(Lmax) + age*(I(Lmax-1));
asymp(Lmax)= (1-sigma)*epsilon*E2(Lmax) + psi*I(Lmax) + FOI(Lmax)*theta*R(Lmax) - (gamma+mu(Lmax))*A(Lmax) + age*(A(Lmax-1))-phi(Lmax)*A(Lmax);
rec(Lmax)=gamma*A(Lmax) - FOI(Lmax)*theta*R(Lmax) - (delta+mu(Lmax))*R(Lmax) + age*(R(Lmax-1))-phi(Lmax)*R(Lmax)+rho*Rv(Lmax);
%vaccinated
susv(Lmax)= kappa*phi(Lmax)*S(Lmax)- (tau*FOI(Lmax)+mu(Lmax))* Sv(Lmax) + delta*Rv(Lmax) + age*(Sv(Lmax-1)) - rho*Sv(Lmax);
exposed1v(Lmax)= kappa*phi(Lmax)*E1(Lmax) + tau*FOI(Lmax)*Sv(Lmax) - (epsilon+mu(Lmax))*E1v(Lmax) + age*(E1v(Lmax-1));
exposed2v(Lmax)= kappa*phi(Lmax)*E2(Lmax) + epsilon*E1v(Lmax) - (epsilon+mu(Lmax))*E2v(Lmax) + age*(E2v(Lmax-1));
infv(Lmax)= sigmaV*epsilon*E2v(Lmax) - (psi+mu(Lmax))*Iv(Lmax) + age*(Iv(Lmax-1));
asympv(Lmax)= kappa*phi(Lmax)*A(Lmax) + (1-sigmaV)*epsilon*E2v(Lmax) + psi*Iv(Lmax) + tau*FOI(Lmax)*theta*Rv(Lmax) - (gamma+mu(Lmax))*Av(Lmax) + age*(Av(Lmax-1));
recv(Lmax)=kappa*phi(Lmax)*R(Lmax) + (1-kappa)*phi(Lmax)*(S(Lmax)+E1(Lmax)+E2(Lmax)+A(Lmax)+R(Lmax))+ gamma*Av(Lmax) - tau*FOI(Lmax)*theta*Rv(Lmax) - (delta+mu(Lmax))*Rv(Lmax) + age*(Rv(Lmax-1))-rho*Rv(Lmax);


%bookkeeping
Cas(Lmax) = psi*I(Lmax)+psi*Iv(Lmax);
AllCases(Lmax)=epsilon*(E2(Lmax)+E2v(Lmax));

%bookkeeping for vaccination
VAC=phi.*(S+E1+E2+A+R ).';

Equations=[mat.';sus.';exposed1.';exposed2.';inf.';asymp.';rec.';susv.';exposed1v.';exposed2v.';infv.';asympv.';recv.';Cas.';VAC.';AllCases.'];
end

