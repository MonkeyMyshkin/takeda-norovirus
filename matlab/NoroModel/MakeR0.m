function [ R0 ] = MakeR0(Param,PopulationSize,ContactMatrix, mu, ageGroupBreaks)
%%MakeR0 takes current parameter values and calculates R0

%INPUTS
%Param=vector of alpha, q, omega1, omega2, nu, delta, epsilon, sigma,
%       psi, gamma
%PopulationSize=vector of agegroup sizes
%ContactMatrix=matrix of contact probabilities
%mu=vector of death rates per age
%ageGroupBreaks=vector of age group divisions

%OUTPUTS
%R0

%%
Lmax=length(mu);
%Define parameter names
alpha=Param(1);
q=Param(2);
nu=Param(5);
delta=Param(6);
epsilon=Param(7);
sigma=Param(8);
psi=Param(9);
gamma=Param(10);

%normalise population size
PopulationSize=PopulationSize/(sum(PopulationSize));
%Transmission matrix
%mostly zeros
Transmission =zeros(4*Lmax);
Transition=zeros(4*Lmax);

for index=1:Lmax
    
    ageGroupIndex=find(histcounts(index,[-1,ageGroupBreaks, Lmax]));
    
    for jndex=1:Lmax
        
        ageGroupJndex=find(histcounts(index,[-1,ageGroupBreaks,Lmax]));
        
        Transmission(index,2*Lmax + jndex) = q*ContactMatrix(ageGroupIndex,ageGroupJndex)  *PopulationSize(jndex); %infected
        
        
        Transmission(index, 3*Lmax+ jndex) = q* ContactMatrix(ageGroupIndex,ageGroupJndex) * nu* PopulationSize(jndex); %asymp
        
    end
    
    %Transition matrix
    %lower triangular
    
    
    Transition(index,index)= epsilon+ mu(index);
    Transition(Lmax+index,Lmax+index) = epsilon + mu(index);
    Transition(2*Lmax+index,2*Lmax+index) = psi + mu(index);
    Transition(3*Lmax+index,3*Lmax+index) = gamma + mu(index);
    
    Transition(Lmax+index,index)=-epsilon;
    Transition(2*Lmax+index, Lmax+index)= -sigma*epsilon;
    Transition(3*Lmax+index, Lmax+index) = -(1-sigma)*epsilon;
    Transition(3*Lmax+index, 2*Lmax+index)=-psi;
    
end

NGM=Transmission/Transition;

R0=max(eig(NGM));

end

