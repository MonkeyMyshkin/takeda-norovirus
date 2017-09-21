function [ ProbabilityOfInfection, Vaccinated ,ProbabilityOfInfectionAll] = ProcessCases_Vac( SimulationResult, AgeGroupBreaks )
%PROCESSCASES: takes cumulative case numbers for each week, finds the
%number of cases in each week and stratifies in age groups
Lmax=81;

C=SimulationResult(:,13*Lmax+1:14*Lmax);

Cases(1,:)=C(1,:);

Cases(2:length(SimulationResult(:,1)),:)= max(C(2:end,:)-C(1:end-1,:),0);



%need to find the size of each age group to get the probability

for index=1:Lmax
    ageGroupSize(:,index)=sum(SimulationResult(:,[index,Lmax+index,2*Lmax+index,3*Lmax+index,4*Lmax+index,5*Lmax+index,6*Lmax+index,7*Lmax+index,8*Lmax+index,9*Lmax+index,10*Lmax+index,11*Lmax+index,12*Lmax+index]),2);
end
ProbabilityOfInfectionAll=bsxfun(@rdivide,Cases,sum(ageGroupSize,2));

Cases=AgeStratify(Cases,AgeGroupBreaks);
ageGroupSize1=AgeStratify(ageGroupSize, AgeGroupBreaks);

ProbabilityOfInfection=bsxfun(@rdivide,Cases,sum(ageGroupSize1,2)); %Probability an individual is age a, infected  in week i


%%
%vaccinated part
V=SimulationResult(:,14*Lmax+1:15*Lmax);

Vaccinated(1,:)=V(1,:);

Vaccinated(2:length(SimulationResult(:,1)),:)= max(V(2:end,:)-V(1:end-1,:),0);

Vaccinated=bsxfun(@rdivide,Vaccinated,sum(ageGroupSize1,2));

end

