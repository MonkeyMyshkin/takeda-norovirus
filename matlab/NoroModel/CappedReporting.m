function [ ReportedInfectionNumber] = CappedReporting( ModelOutput, ReportingBaseline, damping )
%%
%CappedReporting implements incidence-proportional reporting

%INPUTS
%ModelOutput= cases per week in each agegroup
%ReportingBaseline=Proportion of cases reported before 2011
%quenching= values for self correcting Markov process: quenching(1)in [0,1]
%           scales quenching(2) for ages up to 37

%OUTPUTS
%ReportedInfectionNumber=Model output with reporting model applied
%%

noAgeGroups=length(ModelOutput(1,:));

for i=1:4     %damping for <37
    ReportedInfectionNumber(1:33+2*52+1,i)=ModelOutput(1:33+2*52+1,i)*ReportingBaseline(i);
    
    ReportedInfectionNumber(33+2*52+1:418,i)=ModelOutput(33+2*52+1:418,i).*ReportingBaseline(i).*exp(-damping(1)*damping(2)*ModelOutput(33+2*52+1:418,i));
end
for i=5:noAgeGroups     %damping for >37
    ReportedInfectionNumber(1:33+2*52+1,i)=ModelOutput(1:33+2*52+1,i)*ReportingBaseline(i);
    
    ReportedInfectionNumber(33+2*52+1:418,i)=ModelOutput(33+2*52+1:418,i).*ReportingBaseline(i).*exp(-damping(2)*ModelOutput(33+2*52+1:418,i));
end


end

