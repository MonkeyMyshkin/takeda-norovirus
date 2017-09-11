function [AcceptReject,AccProp,LogLik] = MCMCstep(Param,omega2,mu,theta,ContactMatrix,x0,AccCurrent,ageGroupBreaks,GermanCaseNotification,PopulationSize,GermanPopulation,Dispersion,k,d,damping)
%%
%MCMCstep Takes proposed particle values and assesses likelihood and Bayes formula

%INPUTS
%param=parameters includig reporting parameters
%omega2=vector of seasonal offsets for each age group
%mu=death rates per age
%theta=susceptibilities of recovered individuals per season
%ContactMatrix=probabilities of contact between age groups
%x0= initial conditions for simulation
%AccCurrent=current value of likehood and prior probabilities
%ageGroupBreaks=vector of age group divisions
%GermanCaseNotification=data from SurvStat
%PopulationSize=population size of each age group
%GermanPopulation=population size of each annual age group
%Dispersion=dispersion parameter for negative binomial likelihood
%k=degree of diagonalisation of contact matrix
%d=degree of assymetry of contact matrix
%damping=terms for selfcorrecting markov process, damping(1) is a scaling
%        term for 0-37 years and damping(2) is a baseline for all years

%OUTPUTS
%AcceptReject=binary:1 ifaccept,0 if reject
%AccProp=Proposedlikelihood+prior probabilties
%LogLik=Log Likelihood value

%%
%Simulate each season
[ ~,SimulationResult ] =SimulateSeasons(Param(1:9),omega2,mu,theta,ContactMatrix,x0);

%stratify data
[ StratifiedCases ] = AgeStratify( GermanCaseNotification, ageGroupBreaks );

%stratify model output into age groups
[ ModelOutput ] =  ProcessCases( SimulationResult, ageGroupBreaks );        %stratify infected individuals
ModelOutput=ModelOutput*sum(PopulationSize);

%REPORTING MODEL
[ ReportingBaseline ] = Param(10)*[1 Param(11)*ones(1,2) Param(12) Param(13) Param(14) Param(15)]; 
[ ReportedInfectionNumber] = CappedReporting( ModelOutput, ReportingBaseline, damping );

%LIKELIHOOD
[ LogLik ] = NBLikelihood( StratifiedCases, ReportedInfectionNumber, Dispersion);

%PRIOR PROBABILITIES
PriorProbs=Priors(Param(1:9),theta,Param(10:15),Dispersion,k,d,damping,GermanPopulation,ContactMatrix,mu,ageGroupBreaks);

%ACCEPTANCE
AccProp=LogLik+PriorProbs;
[AcceptReject] = ACCEPT( AccProp, AccCurrent);

end

