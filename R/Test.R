source('PrepareData.R')

#ODE Parameters to estimate
params = c(1.0/(0.5*365),1.5,0.1,0.1,1/(10*365),1.0/(30*6))
names(params) <- c('alpha','q','omega','nu','delta','gamma')

# Seasonal ajdustment (set to 0)
theta = rep(0,9)

# Simulate model
x<-SimulateSeasons( c(params[1],params[2],params[3],params[4],params[5],1 ,0.735415 ,0.5,params[6]), B, mu,theta,Cm,Lmax,10)

ProbC <- AgeStratify( x, ageGroupBreaks )
StratifiedSim<-t(t(ProbC)*PopulationSize)

ReportingBaseline = c(1*c(1,1.0,1.0,1.0,1.0,1.0),1.0) 
damping = c(0.0,1e-5)

# Expected cases reports
ReportedInfections  = CappedReporting(StratifiedSim,ReportingBaseline,damping)

# Simulate model and calculate likelihood for test parameters
NBLikelihood(c(params,ReportingBaseline,damping),B,mu,theta,Cm,Lmax,ageGroupBreaks,StratifiedCases)


