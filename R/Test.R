source('PrepareData.R')

params = c(1.0/(0.5*365),1.0,0.1,0.5,1/(10*365),1/2.0,0.5,1/5.0,0.5,0.0)
names(params) <- c('alpha','q','omega','nu','delta','epsilon','sigma','psi','gamma','theta')

theta = 0.0

SimulateSeasons(params, mu,theta,Cm)

#calculate LL and Priors

[ ~,SimulationResult ] =SimulateSeasons(Param(1:9),omega2,mu,theta,ContactMatrix,x0);

%stratify data
[ StratifiedCases ] = AgeStratify( GermanCaseNotification, ageGroupBreaks );

%stratify model output into age groups
[ ModelOutput ] =  ProcessCases( SimulationResult, ageGroupBreaks );        %stratify infected individuals
ModelOutput=ModelOutput*sum(PopulationSize);

%REPORTING MODEL
[ ReportingBaseline ] = [Param(10)*[1 Param(11)*ones(1,2) Param(12) Param(13) Param(14)] Param(15)]; 
[ ReportedInfectionNumber] = CappedReporting( ModelOutput, ReportingBaseline, damping );

