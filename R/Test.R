source('PrepareData.R')
source('InferenceFunctions.R')


#ODE Parameters to estimate
params = c(1.5,0.1,0.1,1/(10*365),1.0/(30*6))
names(params) <- c('q','omega','nu','delta','gamma')

ReportingBaseline = c(1*c(1,1.0,1.0,1.0,1.0,1.0),1.0) 
damping = c(0.0,1e-5)


# Seasonal ajdustment (set to 0)
theta = rep(0,9)

# Simulate model
x<-SimulateSeasons( c(exp(-5.98405),params[1],params[2],params[3],params[4],1 ,0.735415 ,0.5,params[5]), B, mu,theta,Cm,Lmax,10)

ProbC <- AgeStratify( x, ageGroupBreaks )
StratifiedSim<-t(t(ProbC)*PopulationSize)



# Expected cases reports
ReportedInfections  = CappedReporting(StratifiedSim,ReportingBaseline,damping)

# Simulate model and calculate likelihood for test parameters
NBLikelihood(c(params,ReportingBaseline,damping),B,mu,theta,Cm,Lmax,ageGroupBreaks,StratifiedCases)


# Particle Filter
#specifics
noSeasons = 9    #number of seasons
noParam = 14     #number of parameters to be estimated
noParticles=100   #number of particles

# Initial covariance matrix for multivariate normal proposal distribution
init_proposal=(0.1^2)*diag(rep(1,noParam))/noParam;

require('Matrix')
require('parallel')
require('mnormt')

init_pop <- lapply(1:noParticles,function(x){sample_prior()})
out_list<-mclapply(init_pop,NBLikelihood,mc.cores=4,B,mu,theta,Cm,Lmax,ageGroupBreaks,StratifiedCases)
particle_likelihood <- unlist(out_list)
particle_weights <- rep(1,noParticles)

out_list<-filter_particles(init_pop, particle_likelihood, particle_weights, init_proposal,adapt=FALSE)
