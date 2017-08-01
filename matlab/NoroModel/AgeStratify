function [ StratifiedCases ] = AgeStratify( Cases, ageGroupBreaks )
%%
%AgeStratify groups cases into by age group from annual counts
%Inputs: Cases= a vector or matrix of case numbers for each annual age
%               group
%        ageGroupBreaks= a vector of age group divisions 
%Outputs: StratifiedCases=counts grouped by age group
%%
maxAge=length(Cases(1,:));

%need upper bound on ages as limit for collecting counts
ageGroupBreaks=[ageGroupBreaks, maxAge];

%Calculate the number of age groups
noAgeGroups=length(ageGroupBreaks);

%for each week, group cases into each age group ignoring NaN
StratifiedCases(:,1)= nansum(Cases(:, 1 : ageGroupBreaks(1)),2);
for index = 2 : noAgeGroups
    StratifiedCases(:,index)= nansum(Cases(:, ageGroupBreaks(index-1)+1 : ageGroupBreaks(index)),2);
end
    
end

