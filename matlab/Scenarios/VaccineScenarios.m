addpath(genpath('C:/Users/kamg2/Documents/takeda-norovirus'))
%Program to random sample vaccine scenarios and measure: cases averted,
%outpatient visits averted, hospitalisations and deaths averted per dose
%per year
Lmax=81;

%scenario space is 8 dimensional:
% Efficacy, Coverage, Age groups targeted, Vaccine duration, Symptom
% reduction, Infectiousness reduction, Susceptibility reduction, immunizing
Efficacy=[0.5 0.7 0.9];

Coverage=[0.5 0.7 0.9]/365; %coverage of the population in a year

A0=[1;zeros(Lmax-1,1)];
A1=[0;1;zeros(Lmax-2,1)];
A2=[0;0;1;zeros(Lmax-3,1)];
A3=[zeros(50,1); ones(20,1);zeros(Lmax-70,1)];
A4=[zeros(70,1); ones(11,1)];
A5=zeros(Lmax,1); A5([66,71,76,81])=1;
A6=[zeros(50,1); ones(31,1)];  
AgeGroupsTargeted=[A0, A1, A2, A3, A4, A5, A6];

VaccineDuration=[1/365, 1/(3*365), 1/(5*365)];

SymptomReduction=[0 0.3 0.5 0.7 1]; %this scales the proportion symptomatic

InfectiousnessReduction=[0 0.3 0.5 0.7 1];

SusceptibilityReduction=[0 0.3 0.5 0.7 1];

Sterilizing=[1 0];

%proportion outpatient, hospitilized or dead
outpatient=[0.168*ones(1,5), 0.168*ones(1,10), 0.062*ones(1,10), 0.064*ones(1,20), 0.054*ones(1,20), 0.103*ones(1,16)];
hospitilized=[0.00428*ones(1,5), 0.00182*ones(1,13), 0.00228*ones(1,47), 0.01733*ones(1,16)];
fatalities=[6.25e-6*ones(1,5),4.66e-6*ones(1,60),4.35e-4*ones(1,16)];

%For a number of samples of scenario space, calculate reduction in cases
NoScenarios=100;
NoSamples=1;

tic
parfor index=1:NoScenarios
    AgeGrp=randsample(7,1);    Eff=randsample(Efficacy,1);     Cov=randsample(Coverage,1);
    phi=Eff*Cov*AgeGroupsTargeted(:,AgeGrp);
    kappa=randsample(Sterilizing,1);
    tau=randsample(SusceptibilityReduction,1);
    zeta=randsample(InfectiousnessReduction,1);
    SigmaV=randsample(SymptomReduction,1);
    rho=randsample(VaccineDuration,1);
    
    %sampleIndices=randsample(length(PARTICLE(:,1)),NoSamples);
    %(PARTICLE(sampleIndices(jndex),:))
    for jndex=1:NoSamples
        
        [Red,Red_asymp,Vac]=ReductionInCases( PopulationSize,ageGroupBreaks, ContactAgesPerAgeGroup, ...
            T, median(PARTICLE),omega2, mu, phi.',rho, kappa, tau,zeta,SigmaV,0 );
        
        Reduction(index,jndex)=sum(Red);
        Reduction_asymp(index,jndex)=sum(Red_asymp);
        %RepReduction(index,jndex)=RepRed;
        Vaccinated(index,jndex)=sum(Vac(:));
        
        %outpatient, hospitilization and mortality \cite{Bartsch2012}
        OUTPATIENT(index,jndex)=sum(Red.*outpatient);
        HOSPITILIZED(index,jndex)=sum(Red.*hospitilized);
        FATALITIES(index,jndex)=sum(Red.*fatalities);
    end
    
    Scenario(index,:)=[AgeGrp, Eff, Cov, kappa, tau, zeta, SigmaV, rho];
    
end
Reduction=Reduction(:); Reduction_asymp=Reduction_asymp(:); Vaccinated=Vaccinated(:);
OUTPATIENT=OUTPATIENT(:);   HOSPITILIZED=HOSPITILIZED(:);   FATALITIES=FATALITIES(:);
Scenario=repmat(Scenario,NoSamples,1);

time=toc
save('VACCINESCENARIOS')
%% symptomatic

AUC=PlotVaccinationScenarios( Scenario,Reduction./Vaccinated);
title('Symptomatic')
%% outpatient

AUC_out=PlotVaccinationScenarios( Scenario,OUTPATIENT./Vaccinated);
title('Outpatients')
%% hospitilized
AUC_hosp=PlotVaccinationScenarios( Scenario,HOSPITILIZED./Vaccinated);
title('Hospitilization')
%% dead

AUC_dead=PlotVaccinationScenarios( Scenario,FATALITIES./Vaccinated);
title('Deaths')
%%
AUC_asymp=PlotVaccinationScenarios( Scenario,Reduction_asymp./Vaccinated);
h=get(gcf,'children');
set(h,'FontSize',16)
title('All cases')
%%
indices=find(Scenario(:,1)==1);
AUC=PlotVaccinationScenarios( Scenario(indices,:),Reduction_asymp(indices)./Vaccinated(indices));
h=get(gcf,'children');
set(h,'FontSize',16)
%%