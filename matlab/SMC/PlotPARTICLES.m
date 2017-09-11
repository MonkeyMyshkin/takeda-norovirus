function [NoisyModel] = PlotPARTICLES( index, x, z,justmean)
%function to plot output for certain chain inindex

%define parameters at this point
ParamFinal=(x(index,1:15));
dispFinal=(x(index,16));

%susceptibility profile
theta=z(index,:);

%contact parameters
kFinal=(x(index,17));
dFinal=(x(index,18));

%damping
damp=(x(index,19:20));

%load data to compare
load('GermanCaseNotification.mat')
load('GermanPopulation.mat')
load('GermanContactData.mat')
load('GermanParticipantData.mat')
listOfIDs=GermanParticipantData(:,1);
%death function, population distribution
%maximum age group in case notifications is 80. Therefore our max age will
%be 80, indexed as 81.
MaxAge=80;
Lmax=MaxAge+1;
noSeasons=9;
%MU
mu=MakeMu(Lmax);    %daily death rate
%AGEGROUPS
ageGroupBreaks=[8 18 26 37 50 70];
noAgeGroups=length(ageGroupBreaks)+1;
PopulationSize=AgeStratify(GermanPopulation.',ageGroupBreaks);
%seasonal offset
omega2=[0.29345  0.37976  0.56963  0.60415  0.82855  1.0724  1.0357];

%generate contact matrix
[ ~,ContactAgesPerAgeGroup,T] = MakeContactMatrix( ageGroupBreaks, listOfIDs , GermanContactData , GermanParticipantData );
[ NewContactMatrix ] = ContactTwist( ContactAgesPerAgeGroup, T, kFinal, dFinal );

%initial conditions
x0=MakeInitialConditions( ParamFinal,omega2,theta,mu , NewContactMatrix );

%Simulate epidemic
[ ~,SimulationResult ] =SimulateSeasons(ParamFinal(1:9),omega2,mu ,theta, NewContactMatrix , x0);

%stratify data
[ StratifiedCases ] = AgeStratify( GermanCaseNotification, ageGroupBreaks );

%stratify model output like data
[ ModelOutput ,~] =  ProcessCases( SimulationResult, ageGroupBreaks );        %stratify infected individuals
ModelOutput=ModelOutput.*sum(PopulationSize);

%ReportingParam is magnitude of reporting factor
[ ReportingBaseline ] = ParamFinal(10)*[1 ParamFinal(11)*ones(1,2) ParamFinal(12) ParamFinal(13) ParamFinal(14) ParamFinal(15)];
[ ReportedInfectionNumber] = CappedReporting( ModelOutput, ReportingBaseline, damp );

noSamples=1;
%generate negative binomial noise
%with overdispersion %define reported infection number as mean
P=repmat(bsxfun(@rdivide,bsxfun(@times,dispFinal,ReportedInfectionNumber),bsxfun(@plus,1,bsxfun(@times,dispFinal,ReportedInfectionNumber)) ),[1,1,noSamples]); %probability of fail

R=1/dispFinal; 
NoisyModel=nbinrnd(R,1-P);

%plot cases
for index=1:length(ageGroupBreaks)+1
    subplot(4,2,index)
    hold on
    if justmean==1
        ModelOutput=ReportedInfectionNumber(:,index)/PopulationSize(index);
        p=plot(ModelOutput,'LineWidth',1,'color','m');
        p.Color(4)=0.5;
    else
        for jndex=1:noSamples
            p=plot(NoisyModel(:,index,jndex).'/PopulationSize(index),'b');
            p.Color(4)=0.1;
        end
    end
end

for index=1:length(ageGroupBreaks)+1
    subplot(4,2,index)
    hold on
    plot(StratifiedCases(:,index)/PopulationSize(index),'LineWidth',2,'color','k');
end
%% serology
% figure
% %Seropositive
% for index=1:Lmax
%     ageGroupSize(:,index)=sum(SimulationResult(:,[index,Lmax+index,2*Lmax+index,3*Lmax+index,4*Lmax+index,5*Lmax+index,6*Lmax+index]),2);
% end
% for i=1:Lmax
%     seropositive(i)=1-mean(SimulationResult(:,Lmax+i))/mean(ageGroupSize(:,i));
% end
%
% load('serologicalPelosi.mat')
% load('serologicalSon.mat')
% serology=[serologicalPelosi;serologicalSon];
% scatter(serology(:,1),serology(:,3)./serology(:,2),10*serology(:,2))
% hold on
% plot(seropositive)
end

