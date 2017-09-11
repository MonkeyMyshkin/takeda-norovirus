function [Equations ] = MakeMSEIRSequations(q,C, mu, nu, gamma , delta, alpha)
%Makes equations for simple SIRS model with age structure
%%
%INPUTS
%q          transmission probability given contact
%C          contact matrix
%mu         age dependent death rate
%nu         rate of loss of immunity
%gamma      recovery rate
%delta      rate of loss of maternal antibodies
%alpha      rate individuals move from latent to infectious

%OUTPUTS
%Equations  System of symbolic model equations

%%


%number of age groups
Lmax=length(mu);

%symbolic variables
M=sym('M',[Lmax,1]);S=sym('S',[Lmax,1]);E=sym('E',[Lmax,1]);I=sym('I',[Lmax,1]);R=sym('R',[Lmax,1]);

%force of infection for each age group
for i=1:Lmax
    FOI(i)=sum(q*C(i,:)*I);
end

%age group 1
mat(1)=-(delta +mu(1))*M(1);                        %Maternal antibody class
sus(1)=delta*M(1)-(FOI(1)+mu(1))* S(1) + nu*R(1);   %Susceptible
exposed(1)=FOI(1)*S(1) -(alpha+mu(1))*E(1);         %Exposed or latent
inf(1)=alpha*E(1)-(gamma+mu(1))*I(1);               %Infected
rec(1)=gamma*I(1)-(nu+mu(1))*R(1);                  %Recovered
dead(1)=mu(1)*(M(1)+S(1)+E(1)+I(1)+R(1));           %Dead

%age groups 2+
for i=2:Lmax
    mat(i)=-(delta +mu(i))*M(i);
    
    sus(i)=delta*M(i)-(FOI(i)+mu(i))* S(i) + nu*R(i);
    
    exposed(i)= FOI(i)*S(i)-(alpha+mu(i))*E(i);
    
    inf(i)= alpha*E(i)- (gamma+mu(i))*I(i);
    
    rec(i)=gamma*I(i)-(nu+mu(i))*R(i);
    
    dead(i)=mu(i)*(M(i)+I(i)+S(i)+E(i)+R(i));
end

%concatenate equations
Equations=[mat.';sus.';exposed.';inf.';rec.';dead.'];
end

