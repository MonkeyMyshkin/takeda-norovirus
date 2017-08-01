MATLAB CODE: MAP OF FUNCTIONS AND FUNCTION DESCRIPTIONS

MAP OF FUNCTIONS

NOROMODELSMC
•	AgeStratify
•	MakeContactMatrix
•	ContactTwist
•	MCMCstep
	o	SimulateSeasons*
	o	AgeStratify
	o	ProcessCases
			AgeStratify
	o	CappedReporting
	o	NBLikelihood
	o	Priors
			MakeR0
	o	ACCEPT
•	ESS
•	validateCovMatrix
•	MakeIntialConditions

*SimulateSeasons depends on ODE model function created using script

BUILDMODELODE
	o	MakeEquations
	
	
FUNCTION DESCRIPTIONS

AGESTRATIFY
INPUTS
	Cases: 						A vector or matrix of counts for each annual age 
	ageGroupBreaks: 			A vector of age group divisions
OUTPUTS
	StratifiedCases: 			A vector of grouped cases for each age group in ageGroupBreaks

MAKECONTACTMATRIX
INPUTS
	ageGroupBreaks
	GermanContactData:			A matrix from POLYMOD data where first column is a vector of GlobalIDs and the second column is a vector of contact ages
	GermanParticipantData 		A matrix from POLYMOD data where first column is a vector of GlobalIDs, the second is a vector of PartipantAges and the third indicates whether the survey was completed on a weekday
OUTPUTS
	contactMatrix 				A symmetrised average probability of contacts between each age group
	ContactAgesPerAgeGroup		A raw matrix of contacts between each age group
	T 							Total numbers of participants in each age group

CONTACTTWIST
INPUTS
	ContactAgesPerAgeGroup
	T
	k							Degree of diagonalization of contact matrix, k=0 suggests C=Identity, k=1 suggests C is unchanged
	d 							Degree of asymmetry of contact matrix, d=0 suggests matrix is symmetric, d=1 suggests matrix is asymmetric
OUTPUTS
	TwistedContactMatrix		Partially symmetrised, partially diagonalised contact matrix

MAKEEQUATIONS
INPUTS
	Param						Parameter vector including alpha, q, omega1, omega2, nu, delta, epsilon, sigma, psi, gamma, theta and chi. Symbolic.
	mu							Vector of death rates per annual age group. Symbolic.
	theta						Susceptibility of recovered individuals. Symbolic scalar at this stage.
	StratifiedContactMatrix 	Matrix of contact probabilities between age groups. Symbolic.
	ageGroupBreaks
OUTPUTS
	Equations					Differential equations for each age and disease state.

SIMULATESEASONS
INPUTS
	Param 						Parameters including alpha, q, omega1, omega2, nu, delta, epsilon, sigma, psi and chi. Vector of doubles
	mu
	theta						Vector of per-season susceptibility of recovered individuals. 
	ContactMatrix				Matrix of contact probabilities for each age group.
	x0 							Vector of initial conditions for solving ODEs
OUTPUTS
	TIME						Vector of time points, in weeks, from solver
	SimulationResult			matrix of solutions to ODEs for each age group and week

PROCESSCASES
INPUTS
	SimulationResult
	ageGroupBreaks
OUTPUTS
	ProbabilityOfInfection		Proportion of each age group symptomatically infected in each week

CAPPEDREPORTING
INPUTS
	ModelOutput					Cases per week in each age group. Calculated in MCMCstep as ProbabilityOfInfection*ageGroupSize
	ReportingBaseline			Flat reporting rate before 2011. Vector of doubles.
	Damping						Self-correcting Markov Process parameters. Damping (1) in [0,1], is scaling factor of damping (2) for age groups under 37.  
								This ensures that older agegroups experience more damping than younger.
OUTPUTS
	ReportedInfectioonNumber	Model output with reporting model applied

NBLIKELIHOOD
INPUTS
	StratifiedCases				Observed cases from Survstat stratified into same age groups as model output. Matrix of doubles.
	ReportedInfectionNumber
	Dispersion					Dispersion parameter for negative binomial likelihood. Scalar double.
OUTPUTS
	LL							Log likelihood

PRIORS
INPUTS
	Param						Parameters including alpha, q, omega1, omega2, nu, delta, epsilon, sigma, psi and chi. Vector of doubles
	theta						Vector of per-season susceptibility of recovered individuals. 
	ReportingParam				Vector of baseline reporting values. Value is shared for 18-26 and 26-37 year olds
	Dispersion					Dispersion parameter for negative binomial likelihood
	k							Degree of diagonalization of contact matrix
	d							Degree of asymmetry in contact matrix
	damping						Self-correcting Markov Process parameters. Damping (1) in [0,1], is scaling factor of damping (2) for age groups under 37.  
								This ensures that older agegroups experience more damping than younger.
	GermanPopulation			Vector of age group sizes for German Population
	ContactMatrix				Matrix of contact probabilities for given parameters.
	mu
	ageGroupBreaks
OUTPUTS
	PriorProbabilities			Sum of log prior probabilities for current parameters.

MAKER0
INPUTS
	Param
	PopulationSize
	ContactMatrix
	mu
	ageGroupBreaks
OUTPUTS
	R0

ACCEPT
INPUTS
	accProp						Proposed likelihood and prior probability
	accCurrent					Current likelihood and prior probability
OUTPUTS
	accept						Binary term-1 if accept and 0 if reject

ESS
INPUTS
	weight						Vector of particle weights before resampling
OUTPUTS
	ESS							Effective sample size

MAKEINITIALCONDITIONS
INPUTS
	Param
	theta
	mu
	ContactMatrix
OUTPUTS
	InitialConditions			Initial conditions,x0, for SimulateSeasons

MCMCSTEP
INPUTS
	param
	mu
	theta
	ContactMatrix
	x0
	AccCurrent
	ageGRoupBreaks
	GermanCaseNotification		SurvStat data
	PopulationSize
	GermanPopulation
	Dispersion
	k
	d
	damping
OUTPUTS
	AcceptReject				Binary output, 1 for accept, 0 for reject
	AccProp						Likelihood+ prior probabilities
	LogLik						Log Likelihood for proposed

MAKEMU
INPUTS
	Lmax 						Maximum age index
OUTPUTS
	mu							Death rate per age
