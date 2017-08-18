function [ ProbabilityOfInfection ] = ProcessCases( SimulationResult, AgeGroupBreaks )
%%
%PROCESSCASES: takes cumulative case numbers for each week, finds the
%number of cases in each week and stratifies in age groups

%INPUTS
%SimulationResult=solution of ODEs for each week in TIME and each age
%ageGroupBreaks=a vector of age group divisions 

%OUTPUTS
%ProbabilityOfInfection=proportion of age group that is infected in each
%                       week
%%
Lmax=81;

%cumulative numbers of symptomatic infections per week
C=SimulationResult(:,8*Lmax+1:9*Lmax);

%cases that occurred in that week are the difference betwee weeks
Cases(1,:)=C(1,:);
Cases(2:length(SimulationResult(:,1)),:)= max(C(2:end,:)-C(1:end-1,:),0);

%age stratify cases
Cases=AgeStratify(Cases,AgeGroupBreaks);

%calculate the size of each age group to get the probability
for index=1:Lmax
    ageGroupSize(:,index)=sum(SimulationResult(:,[index,Lmax+index,2*Lmax+index,3*Lmax+index,4*Lmax+index,5*Lmax+index,6*Lmax+index,7*Lmax+index]),2);
end
ageGroupSize=AgeStratify(ageGroupSize, AgeGroupBreaks);

ProbabilityOfInfection=bsxfun(@rdivide,Cases,sum(ageGroupSize,2)); %Probability an individual is age a, infected  in week i

end

