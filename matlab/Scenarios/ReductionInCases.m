function [ Reduction,Reduction_asymp , Vacc] = ReductionInCases(PopulationSize,ageGroupBreaks, ContactAgesPerAgeGroup, T, PARTICLE,omega2, mu, phi, rho, kappa, tau,zeta, sigmaV, figDisplay )
%function to calculate % reduction in cases over time period for parameter
%set
Lmax=81;

%DEFINE PARAMETERS
Param=PARTICLE(1:9);   %alpha, q, omega1,  nu, delta, epsilon, sigma, psi, gamma, chi
ReportingBaseline = [PARTICLE(10) PARTICLE(11)*ones(1,2) PARTICLE(12) PARTICLE(13) PARTICLE(14) PARTICLE(15)];
k=PARTICLE(17);
d=PARTICLE(18);
quenching=PARTICLE(19:20);
theta=PARTICLE(21:end);

%contact and ICs
[ NewContactMatrix ] = ContactTwist( ContactAgesPerAgeGroup, T, k, d );
x0=MakeInitialConditions( Param,omega2,theta,mu , NewContactMatrix );

%SIMULATE WITHOUT VACCINATION
[ ~,SimulationResult ] =SimulateSeasons(Param(1:9),omega2 ,mu, theta, NewContactMatrix , x0);
[ ModelOutput ,ModelOutputAll ] =  ProcessCases( SimulationResult, ageGroupBreaks );
[ ModelOutput_asymp ] =  ProcessCases_asymp( SimulationResult, ageGroupBreaks );

%Vaccine simulation and parameters
sigmaV=sigmaV*Param(7);   %proportion of infected vaccianted individuals who are symptomatic
ParamVac=[rho, kappa,tau, zeta, sigmaV];

%SIMULATE WITH VACCINATION
[~,SimulationResultVac ] = SimulateSeasons_Vaccination( Param,ParamVac,omega2,mu ,phi,theta, NewContactMatrix , x0);
%[~,SimulationResultVac ] = SimulateLifetime_Vaccination( Param,ParamVac,mu ,phi,mean(theta), NewContactMatrix , x0);
[ ModelOutputVac,Vacc ,ModelOutputVacAll] =  ProcessCases_Vac( SimulationResultVac, ageGroupBreaks );
[ ModelOutputVac_asymp,~] =  ProcessCases_Vac_asymp( SimulationResultVac, ageGroupBreaks );

%SCALE UP TO POPULATION LEVEL
PopModelOutput=ModelOutput*sum(PopulationSize);
PopModelOutputAll=ModelOutputAll*sum(PopulationSize);
PopModelOutputVacAll=ModelOutputVacAll*sum(PopulationSize);
PopModelOutputVac=ModelOutputVac*sum(PopulationSize);
PopModelOutput_asymp=ModelOutput_asymp*sum(PopulationSize);
PopModelOutputVac_asymp=ModelOutputVac_asymp*sum(PopulationSize);

load('GermanPopulation.mat')
Vacc=Vacc.*sum(PopulationSize)/418 *52;

%OPTIONAL FIGURE
if figDisplay
    for index=1:7
        subplot(4,2,index)
        plot(PopModelOutputAll(:,index),'.')
        hold on
        plot(PopModelOutputVac(:,index))
    end
end

%NUMBER OF CASES AVERTED
Reduction=sum(PopModelOutputAll-PopModelOutputVacAll)/418 *52;
Reduction_asymp=sum(PopModelOutput_asymp-PopModelOutputVac_asymp)/418 *52;

%REPORTED CASES
[ ReportedInfectionNumber] = CappedReporting( PopModelOutputAll, ReportingBaseline, quenching );
[ ReportedInfectionNumberVac] = CappedReporting( PopModelOutputVac, ReportingBaseline, quenching );



end

