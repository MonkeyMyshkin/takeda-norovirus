source('truncated.R')

CappedReporting <- function(ModelOutput, ReportingBaseline, damping)
{

# CappedReporting implements incidence-proportional reporting

# INPUTS
# ModelOutput= cases per week in each agegroup
# ReportingBaseline=Proportion of cases reported before 2011
# quenching= values for self correcting Markov process: quenching(1)in [0,1]
# scales quenching(2) for ages up to 37

# OUTPUTS
# ReportedInfectionNumber=Model output with reporting model applied

noAgeGroups = dim(ModelOutput)[2]

ReportedInfectionNumber = ModelOutput*0

for( i in 1:4) #damping for <37
{
 ReportedInfectionNumber[1:(33+2*52+1),i] = ModelOutput[1:(33+2*52+1),i] * ReportingBaseline[i]
 ReportedInfectionNumber[(33+2*52+1):418,i] = ModelOutput[(33+2*52+1):418,i] * ReportingBaseline[i]*exp(-damping[1]*damping[2]*ModelOutput[(33+2*52+1):418,i])
}

for (i in c(5:noAgeGroups))    #damping for >37
{
    ReportedInfectionNumber[1:(33+2*52+1),i]=ModelOutput[1:(33+2*52+1),i]*ReportingBaseline[i]
    ReportedInfectionNumber[(33+2*52+1):418,i]=ModelOutput[(33+2*52+1):418,i]*ReportingBaseline[i]*exp(-damping[2]*ModelOutput[(33+2*52+1):418,i]);
}

return(ReportedInfectionNumber)

}



NBLikelihood <- function(params,B,mu,theta,ContactMatrix,MaxAge,ageGroupBreaks,StratifiedCases)
{

p_sim = c(exp(-5.98405),params[1],params[2],params[3],params[4],1 ,0.735415 ,0.5,params[5])
#names(p_sim) <- c('alpha','q','omega','nu','delta','epsilon','sigma','psi','gamma')

x<-SimulateSeasons(p_sim,B, mu,theta,Cm,Lmax,10)

ProbC <- AgeStratify( x, ageGroupBreaks )
StratifiedSim<-t(t(ProbC)*PopulationSize)

ReportingBaseline = c(params[6]*c(1,rep(params[7],2),params[8],params[9],params[10]),params[11]) 

Dispersion = params[12]

damping = c(params[13],params[14])

ReportedInfectionNumber  = CappedReporting(StratifiedSim,ReportingBaseline,damping)

sum(dnbinom(x=StratifiedCases,mu=ReportedInfectionNumber,size=Dispersion,log=TRUE))

}

sample_prior <- function()
{

 prob = numeric(14)
 
 # q
 prob[1] = runif(1,50,400)
 # omega
 prob[2] = rtrunc(1,'gamma',shape=0.15,scale=1,a=0.01,b=0.3)
 # nu
 prob[3] = rtrunc(1,'gamma',shape=1,scale=1,a=0,b=1)
 # delta
 prob[4] = rtrunc(1,'gamma',shape=1/(2*365),scale=1,a=1/(10*365),b=2/365)
 # gamma
 prob[5] = rtrunc(1,'gamma',shape=1,scale=1,a=0,b=1)
 # Reporting
 prob[6:11] = runif(6,0,1)
 # Dispersion
 prob[12] = rtrunc(1,'gamma',shape=1,scale=1,a=0,b=0.5)
 # Damping
 prob[13] = runif(1,0,1)
 prob[14] = rgamma(1,shape=1e-1,scale=1e-2)
 
 return(prob) 

}

prior_prob <- function(params)
{

prob = numeric(14)

 # q
 prob[1] = dunif(params[1],50,400,log=TRUE)
 # omega
 prob[2] = dtrunc(params[2],'gamma',shape=0.15,scale=1,a=0.01,b=0.3,log=TRUE)
 # nu
 prob[3] = dtrunc(params[3],'gamma',shape=1,scale=1,a=0,b=1,log=TRUE)
 # delta
 prob[4] = dtrunc(params[4],'gamma',shape=1/(2*365),scale=1,a=1/(10*365),b=2/365,log=TRUE)
 # gamma
 prob[5] = dtrunc(params[5],'gamma',shape=1,scale=1,a=0,b=1,log=TRUE)
 # Reporting
 prob[6:11] = dunif(params[6:11],0,1,log=TRUE)
 # Dispersion
 prob[12] = dtrunc(params[12],'gamma',shape=1,scale=1,a=0,b=0.5,log=TRUE)
 # Damping
 prob[13] = dunif(params[13],0,1,log=TRUE)
 prob[14] = dgamma(params[14],shape=1e-1,scale=1e-2,log=TRUE)
 
 return(prob) 
}


