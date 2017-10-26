addpath(genpath('C:/Users/kamg2/Documents/takeda-norovirus'))
%script to run particle filter

%%%% Preliminaries %%%%
%clear workspace
clear all

%DATA
load('GermanCaseNotification.mat')  %survstat data
load('GermanPopulation.mat')        %German population sizes per age
load('GermanContactData.mat')       %Contact data from POLYMOD
load('GermanParticipantData.mat')   %participant data from POLYMOD

%AGE GROUPS
ageGroupBreaks=[8 18 26 37 50 70];
noAgeGroups=length(ageGroupBreaks)+1;
PopulationSize=AgeStratify(GermanPopulation.',ageGroupBreaks); %calcul;ate population sizes for these age groups

%CONTACT MATRIX
[ contactMatrix,ContactAgesPerAgeGroup,T] = MakeContactMatrix( ageGroupBreaks, GermanContactData , GermanParticipantData );

%DEATH FUNCTION
%maximum age group in case notifications is 80. Therefore our max age will
%be 80, indexed as 81.
MaxAge=80;
Lmax=MaxAge+1;
mu=MakeMu(Lmax);    %daily death rate

%randomise start point for any pseudo random number generator
rng('shuffle')

%specifics
noSeasons=9;    %number of seasons
noParam=14;     %number of parameters to be estimated
NumberParticles=1000;

%seasonal offset for each age group
load('omega2.mat')

%PROPOSAL DISTRIBUTION
proposal=(0.1^2)*eye(noParam)/noParam;

%%%% Generate Particles from priors %%%%
parfor particleIndex=1:NumberParticles
    %start conditions for each particle
    alphastart=exp(-5.98405);   %lognrnd(-5.98405,0.185685);
    qstart=(unifrnd(50,400));
    %omega1
    pd=makedist('gamma','a',0.15,'b',1)
    pd=truncate(pd,0.01,0.3);  %from prior
    omega1start=random(pd);
    %nu
    pd=makedist('gamma','a',1,'b',1)
    pd=truncate(pd,0,1);
    nustart=random(pd);    %from prior
    %delta
    pd=makedist('gamma','a',1/(2*365),'b',1)
    pd=truncate(pd,1/(10*365),2/365);
    deltastart=random(pd); %from prior
    
    epsilonstart=1;   %never estimate this
    sigmastart=0.735415;    %FIXED
    psistart=0.5;   %rate symptoms are lost
    %gamma
    pd=makedist('gamma','a',1,'b',1)
    pd=truncate(pd,0,1);
    gammastart=random(pd); %from prior %scaling of symptomatic duration
    ReportingStart=[unifrnd(0,1),unifrnd(0,1),unifrnd(0,1)*ones(1,2), unifrnd(0,1), unifrnd(0,1)];    %reporting 1 is baseline
    
    ParamCurrent=([alphastart, qstart, omega1start, nustart, deltastart, epsilonstart, sigmastart, psistart, gammastart, ReportingStart]);
    
    %dispersion
    pd=makedist('gamma','a',1,'b',1)
    pd=truncate(pd,0,0.5);  %from prior
    DispersionPar=random(pd);   %from prior
    
    %contact structure
    kPar=unifrnd(0.01,2);  %from prior
    dPar=0;    %FIXED
    
    %damping
    dampPar=[unifrnd(0,1),gamrnd(1e-1,1e-2)]; %scaling factor for damping parameters between 0 and 1 because damping is less for the young
    
    %susceptibility of recovered individuals
    thetaPar = zeros(1,noSeasons);  %NONE, FIXED
    
    %form together parameters for PARTICLE
    PARTICLE(particleIndex,:)=[ParamCurrent DispersionPar kPar dPar dampPar thetaPar];
    
end

%this is MATLAB quirk, closing the parallel pool can help with potential
%memory leaks
delete(gcp)
reset(symengine)


%%% define equal weights %%%%
Weights=ones(1,NumberParticles);
IterIndex=1;

tic
%%%% Start filter loop %%%%
for IterIndex=IterIndex:100
    
    PARTICLEhistory(:,:,IterIndex)=PARTICLE;
    
    %%%% Normalize weights %%%%
    Weights=Weights/nansum(Weights);
    Weights(isnan(Weights))=0;
    
    
    %%%% RESAMPLING %%%%
    
    %%%% Resample particles using weights %%%%
    if ESS(Weights)/NumberParticles < 0.7
        %%%% Resample particles using weights %%%%
        NewPARTICLEind=datasample(1:NumberParticles,NumberParticles,'Weights',Weights);
        ResampledPARTICLE=PARTICLE(NewPARTICLEind,:);
        Weights=1/NumberParticles*ones(1,NumberParticles);
        resample=1;
    else
        ResampledPARTICLE=PARTICLE;
        resample=0;
    end
    
    
    %%%% Propagate current particles %%%%
    %In order to improve the performane of the parfor loop I slice the
    %particle
    
    %sliced resampled particle
    ParamSliceResampledPARTICLE=ResampledPARTICLE(:,1:15);
    DispSliceResampledPARTICLE=ResampledPARTICLE(:,16);
    kSliceResampledPARTICLE=ResampledPARTICLE(:,17);
    dSliceResampledPARTICLE=ResampledPARTICLE(:,18);
    dampSliceResampledPARTICLE=ResampledPARTICLE(:,19:20);
    thetaSliceResampledPARTICLE=ResampledPARTICLE(:,21:end);
    
    %%%perform MCMC step on each of the particles
    parfor particleIndex=1:NumberParticles
        
        %reset values for each particle
        ParamCurrent=ParamSliceResampledPARTICLE(particleIndex,:);
        DispersionPar=DispSliceResampledPARTICLE(particleIndex);
        kPar=kSliceResampledPARTICLE(particleIndex);
        dPar=dSliceResampledPARTICLE(particleIndex);
        dampPar=dampSliceResampledPARTICLE(particleIndex,:);
        thetaPar=thetaSliceResampledPARTICLE(particleIndex,:);
        
        %%%% Calculate probability of that ResampledParticle for MCMC step %%%%
        
        %contactmatrix and ICS
        [ ContactMatrix ] = ContactTwist( ContactAgesPerAgeGroup, T, kPar ,dPar);
        x0Current=MakeInitialConditions( ParamCurrent,omega2,thetaPar,mu,ContactMatrix);
        
        %calculate LL and Priors
        [ ~, AccCurrent(particleIndex) ,LL(particleIndex)] = MCMCstep(ParamCurrent,omega2,mu,thetaPar,ContactMatrix,x0Current,1,...
            ageGroupBreaks,GermanCaseNotification,PopulationSize,GermanPopulation,DispersionPar,kPar,dPar,dampPar);
        
        %%%% Propose new particle %%%%%
        ParamProp=ParamCurrent; %keep old parameters for fixed parameters
        
        if IterIndex>10 && rand(1)<0.9  %this adapts 90% of the time
            %save the covariance for the transition kernel
            adapt=cov(PARTICLE(:,[2:5,9:12,14:17,19:20]))* 2;
            %indexes correspond to the parameters we are estimating
            
            New=(mvnrnd([ParamCurrent([2:5,9:12,14:15]),DispersionPar,kPar,dampPar],adapt));
            TransitionKernelCovariance=validateCovMatrix(adapt);
            %validateCovMatrix just makes sure no numerical error renders the covariance matrix non-positive-definite
        else
            %otherwise use proposal distribution defined at start
            New=(mvnrnd([ParamCurrent([2:5,9:12,14:15]),DispersionPar,kPar,dampPar],proposal));
            TransitionKernelCovariance=proposal;
        end
        %assign proposed parameters
        ParamProp(2:5)=New(1:4);
        ParamProp(9)=New(5);
        ParamProp(10:12)=New(6:8);
        ParamProp(13)=New(8);
        ParamProp(14)=New(9);
        ParamProp(15)=New(10) ;
        DispersionProp=New(11);
        kProp=New(12);
        dampProp=New(13:14);
        
        %ParamPropParticle
        PRP=[ParamProp DispersionProp kProp dPar dampProp thetaPar];
        
        %%%% Calculate probability of new Particle to compare at accept reject
        [ ContactMatrixProp ] = ContactTwist( ContactAgesPerAgeGroup, T, kProp ,dPar);
        
        %%%% Test proposed values are in range %%%%
        %two conditions- positive values and within prior ranges
        OUTofBOUNDS= min(PRP)<0 || isinf(Priors(ParamProp(1:9),thetaPar,ParamProp(10:15),DispersionProp,...
            kProp,dPar,dampProp,GermanPopulation,ContactMatrixProp,mu,ageGroupBreaks));
        
        if  OUTofBOUNDS  %check positive and in bounds
            AccProp(particleIndex)=-inf;
            LLProp(particleIndex)=-inf;
            AcceptReject=0;
        else
            % Initial conditions
            x0Prop=MakeInitialConditions(ParamProp,omega2,thetaPar,mu,ContactMatrixProp);
            
            %MCMC step
            [ AcceptReject, AccProp(particleIndex) ,LLProp(particleIndex)] = MCMCstep(ParamProp,omega2,mu,thetaPar,ContactMatrixProp,x0Prop,AccCurrent(particleIndex),...
                ageGroupBreaks,GermanCaseNotification,PopulationSize,GermanPopulation,DispersionProp,kProp,dPar,dampProp);
            
        end
        
        
        
        %%%% Update %%%%
        CurrentResampledPARTICLE=[ParamSliceResampledPARTICLE(particleIndex,:),DispSliceResampledPARTICLE(particleIndex),kSliceResampledPARTICLE(particleIndex),...
            dSliceResampledPARTICLE(particleIndex),dampSliceResampledPARTICLE(particleIndex,:),thetaSliceResampledPARTICLE(particleIndex,:)];
        
        %%%% Accept or reject new particle %%%%
        if AcceptReject
            nP= PRP;
            AccCurrent(particleIndex)=AccProp(particleIndex);
            LL(particleIndex)=LLProp(particleIndex);
        else
            nP= CurrentResampledPARTICLE;
            NoAccept(particleIndex)=0;
        end
        
        
        TransitionKernel(particleIndex)=mvnpdf(nP([2,3,4,5,9,10:12,14:15,16,17,19:20]),CurrentResampledPARTICLE([2,3,4,5,9,10:12,14:15,16,17,19:20]),TransitionKernelCovariance);
        newPARTICLE(particleIndex,:)=nP;
    end
    
    PARTICLE=newPARTICLE;
    NumberAccepted(IterIndex)=mean(NoAccept);
    AccCurrent(isinf(AccCurrent))=nan; %ignore infinite values
    %%%% Calculate Weights %%%%
    %IF NOT RESAMPLING AT EVERY STEP
    if resample
        Weights=(AccCurrent+abs(min(AccCurrent))+1)./sum(TransitionKernel);
    else
        Weights= Weights.*((AccCurrent+abs(min(AccCurrent))+1)./sum(TransitionKernel)) ;
    end
    
    Like(IterIndex,:)=LL;
    WeightHistory(IterIndex,:)=Weights/nansum(Weights);
    time(IterIndex)=toc;
    tic
    save('PARTICLEtest')
    delete(gcp)
    reset(symengine)
    fprintf('Iteration is: %d ',IterIndex)
end

