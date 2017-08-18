function [contactMatrix,ContactAgesPerAgeGroup,T] = MakeContactMatrix(ageGroupBreaks,GermanContactData,GermanParticipantData)
%%
%MakeContactMatrix takes Polymod data and produces contact matrix and age
%group counts

%Inputs:
%(POLYMOD data: global ids from contacts and participants, participant and
%contact ages and weekend/weekday info: weekday=1, weekend=0)
%GermanContactData=[GlobalIDs, ContactAges]
%GermanParticipantData=[GlobalIds, ParticipantAges , Weekday?]
%Outputs:
%contactMatrix= a symmetrised average probability of contact between age
%               groups see Bageulin et al. 2013
%ContactAgesPerAgeGroup=a raw matrix of contacts between each age group
%T= Total numbers of participants in each age group

%%
%define a maximum age for age counts
maxAge=99;

ageGroupBreaks=[ageGroupBreaks, maxAge];

noAgeGroups=length(ageGroupBreaks);

noParticipants=length(GermanParticipantData(:,1));

%sorting your participant ages in order
[~,ageSortIndex]=sort(GermanParticipantData(:,2));
GermanParticipantData=GermanParticipantData(ageSortIndex,:);

%for each participant, we note whether they have contact with a contact of
%each yearly age
%DECLARE: row for every participant, column for every age
ContactAgesPerParticipant=zeros(noParticipants,noAgeGroups);

for index= 1 : noParticipants
    %identify all contacts related to current participant
    %isolate the ages of these contacts
    %sort out frequency of contacts in each age group: agegroupbreaks are
    %increased by 1 as they correspond to an index initialised from 1
    ContactAgesPerParticipant(index,:)=histcounts(GermanContactData(GermanContactData(:,1)==GermanParticipantData(index,1),2),[-1,ageGroupBreaks]);
end

%Define the number of participants of each age, T, in each age
%group, W, and number of weekday surveys, Nwd
W=histcounts(GermanParticipantData(:,2),0:maxAge+1);
T=histcounts(GermanParticipantData(:,2),[0,ageGroupBreaks+1]);
Nwd=nnz(GermanParticipantData(:,3));   %number of non zero elements


% Calculate weightings for weekday/weekend participant surveys
% DECLARE: weight diagonal matrix for ease later on
weight=zeros(noParticipants);

for index= 1 : noParticipants
    
    ageGroupIndex= find(ageGroupBreaks>GermanParticipantData(index,2),1);
    
    if GermanParticipantData(index,3)
        weight(index,index) = (W(GermanParticipantData(index,2)+1)/T(ageGroupIndex) )* 5/Nwd;
    else
        weight(index,index) = (W(GermanParticipantData(index,2)+1)/T(ageGroupIndex) )* 2/(noParticipants-Nwd);
    end
    
end

%Update ContactAgesPerParticipant with weighting
ContactAgesPerParticipant=weight*ContactAgesPerParticipant;


%Calculate all contacts for participants in one age group with contacts in
%another
%DECLARE: row , column for every age
ContactAgesPerAgeGroup=zeros(noAgeGroups);

%first need indexes of last participants in each age group
partIndexLast=zeros(1,noAgeGroups);
for index = 1 : noAgeGroups
    partIndexLast(index)=find(GermanParticipantData(:,2)<=ageGroupBreaks(index),1,'last');
end

ContactAgesPerAgeGroup(1,:) = nansum( ContactAgesPerParticipant( 1 : partIndexLast(1), :)) / ...
        trace( weight( 1 : partIndexLast(1), 1 : partIndexLast(1)));


for index = 2 : noAgeGroups
    
    ContactAgesPerAgeGroup(index,:) = sum( ContactAgesPerParticipant( partIndexLast(index-1)+1 : partIndexLast(index), :)) / ...
        trace( weight( partIndexLast(index-1)+1 : partIndexLast(index),partIndexLast(index-1)+1 : partIndexLast(index)));
    
end

%finally symmetrise contactmatrix
%Declare: row, column for every age
contactMatrix=zeros(noAgeGroups);

for index = 1 : noAgeGroups
    for jndex = 1 : noAgeGroups
        
        contactMatrix(index,jndex) = 0.5 * ( ContactAgesPerAgeGroup(index,jndex)/ T(index) + ContactAgesPerAgeGroup(jndex,index)/ T(jndex) );
        
    end
end

end

