addpath(genpath('C:/Users/kamg2/Documents/takeda-norovirus/matlab'))
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
noParam=26;     %number of parameters to be estimated
NumberParticles=50;

%seasonal offset for each age group
omega2=[0.29345  0.37976  0.56963  0.60415  0.82855  1.0724  1.0357];

%PROPOSAL DISTRIBUTION
proposal=(0.1^2)*eye(noParam)/noParam;

%%%% Generate Particles from priors %%%%
parfor particleIndex=1:NumberParticles
    %start conditions for each particle
    alphastart=exp(-5.98405);   %from serology FIXED
    qstart=(unifrnd(10,400));   %arbitrary range- does not have explicit prior
    omega1start=unifrnd(0.05,0.5);  %from prior
    nustart=unifrnd(0,1);           %from prior
    deltastart=normrnd(1/(5.1*365), 5e-5);%unifrnd(1/(30*365), 1);  %broader than informative prior
    epsilonstart=1;   %never estimate this- duration of latency
    sigmastart=unifrnd(0,1); %normrnd(0.735415,0.0960897); %depending on the prior
    psistart=0.5;   %never estimate this- duration of symptoms
    gammastart=gamrnd(1,1); %from priors
    
    ReportingStart=[unifrnd(0,1),unifrnd(0,1),unifrnd(0,1),unifrnd(0,1),unifrnd(0,1), unifrnd(0,1)]; %not full prior in all cases- can adjust
    
    
    ParamCurrent=abs([alphastart, qstart, omega1start, nustart, ...
        deltastart, epsilonstart, sigmastart, psistart, gammastart, ReportingStart]);
    
    %dispersion
    DispersionPar=unifrnd(0,0.4);   %from prior
    
    %contact structure
    kPar=unifrnd(0.01,2);  %from prior
    dPar=unifrnd(0,0.99);    %from prior
    
    %damping
    dampPar=[unifrnd(0,1),unifrnd(9e-7,1e-3)]; %scaling factor for damping parameters between 0 and 1 because damping is less for the young
    
    %susceptibility of recovered individuals
    thetaPar = unifrnd(0,ones(1,noSeasons));
    
    %form together parameters for PARTICLE
    PARTICLE(particleIndex,:)=[ParamCurrent DispersionPar kPar dPar dampPar thetaPar];
    
end

%this is MATLAB quirk, closing the parallel pool can help with potential
%memory leaks
delete(gcp)
reset(symengine)


%%% define equal weights %%%%
Weights=ones(1,NumberParticles);


tic
%%%% Start filter loop %%%%
for IterIndex=1:100
    
    PARTICLEhistory(:,:,IterIndex)=PARTICLE;
    
    %%%% Normalize weights %%%%
    Weights=Weights/nansum(Weights);
    Weights(isnan(Weights))=0;
    
    
    %%%% RESAMPLING %%%%
    
    %%%% Resample particles using weights %%%%
    NewPARTICLEind=datasample(1:NumberParticles,NumberParticles,'Weights',Weights);
    ResampledPARTICLE=PARTICLE(NewPARTICLEind,:);
    Weights=1/NumberParticles*ones(1,NumberParticles);
    
    
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
        
        if IterIndex>1 && rand(1)<0.9  %this adapts 90% of the time
            %save the covariance for the transition kernel
            adapt=cov(PARTICLE(:,[2,3,4,5,7,9,10:15,16,17,18,19:20,21:end]))* 2;
            %indexes correspond to the parameters we are estimating
            
            New=(mvnrnd([ParamCurrent([2,3,4,5,7,9,10:15]),DispersionPar,kPar,dPar,dampPar,thetaPar],adapt));
            TransitionKernelCovariance=validateCovMatrix(adapt);
            %validateCovMatrix just makes sure no numerical error renders the covariance matrix non-positive-definite
        else
            %otherwise use proposal distribution defined at start
            New=(mvnrnd([ParamCurrent([2,3,4,5,7,9,10:15]),DispersionPar,kPar,dPar,dampPar,thetaPar],proposal));
            TransitionKernelCovariance=proposal;
        end
        %assign proposed parameters
        ParamProp(2:5)=New(1:4);
        ParamProp(7)=New(5);
        ParamProp(9)=New(6);
        ParamProp(10:15)=New(7:12);
        DispersionProp=New(13);
        kProp=New(14);
        dProp=New(15);
        dampProp=New(16:17);
        thetaProp=New(18:end);
        
        %ParamPropParticle
        PRP=[ParamProp DispersionProp kProp dProp dampProp thetaProp];
        
        %%%% Calculate probability of new Particle to compare at accept reject
        [ ContactMatrixProp ] = ContactTwist( ContactAgesPerAgeGroup, T, kProp ,dProp);
        
        %%%% Test proposed values are in range %%%%
        %two conditions- positive values and within prior ranges
        OUTofBOUNDS= min(PRP)<0 || isinf(Priors(ParamProp(1:9),thetaProp,ParamProp(10:15),DispersionProp,...
            kProp,dProp,dampProp,GermanPopulation,ContactMatrixProp,mu,ageGroupBreaks));
        
        if  OUTofBOUNDS  %check positive and in bounds
            AccProp(particleIndex)=-inf;
            LLProp(particleIndex)=-inf;
            AcceptReject=0;
        else
            % Initial conditions
            x0Prop=MakeInitialConditions(ParamProp,omega2,thetaProp,mu,ContactMatrixProp);
            
            %MCMC step
            [ AcceptReject, AccProp(particleIndex) ,LLProp(particleIndex)] = MCMCstep(ParamProp,omega2,mu,thetaProp,ContactMatrixProp,x0Prop,AccCurrent(particleIndex),...
                ageGroupBreaks,GermanCaseNotification,PopulationSize,GermanPopulation,DispersionProp,kProp,dProp,dampProp);
            
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
        end
        
        
        TransitionKernel(particleIndex)=mvnpdf(nP([2,3,4,5,7,9,10:15,16,17,18,19:20,21:end]),CurrentResampledPARTICLE([2,3,4,5,7,9,10:15,16,17,18,19:20,21:end]),TransitionKernelCovariance);
        newPARTICLE(particleIndex,:)=nP;
    end
    
    PARTICLE=newPARTICLE;
    
    AccCurrent(isinf(AccCurrent))=nan; %ignore infinite values
    %%%% Calculate Weights %%%%
    Weights=(AccCurrent+abs(min(AccCurrent))+1)./sum(TransitionKernel);
    
    Like(IterIndex,:)=LL;
    WeightHistory(IterIndex,:)=Weights;
    time(IterIndex)=toc;
    tic
    save('PARTICLEtest')
    delete(gcp)
    reset(symengine)
    fprintf('Iteration is: %d ',IterIndex)
end

