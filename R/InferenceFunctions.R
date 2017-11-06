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

ESS <- function(weight)
{
#ESS calculates effective sample size based on weights

# INPUTS
# weight=vector of particle weights 

# OUTPUTS
# ESS=effective samplesize

# Do not assume weights are normalised, so do so here
weight = weight / sum(weight)
return(1/sum(weight^2))
}

filter_particles <- function(particles, particle_likelihood, particle_weights, init_proposal,adapt=FALSE)
{

  # Normalise weights 
  proposed_particles = particles
  particle_weights=particle_weights/sum(particle_weights)
  #weights(isnan(weights))=0

  proposal_cov = init_proposal

  # RESAMPLING 
    
  # If ESS Falls below threshold resample particles using weights 
  if( ESS(particle_weights)/noParticles < 0.7)
  {
        # Resample particles 
        samp_ind = sample(1:noParticles,noParticles,prob=particle_weights,replace=TRUE)
        particles=particles[samp_ind]
        particle_likelihood = particle_likelihood[samp_ind]
        particle_weights=(1/noParticles)*rep(1,noParticles);
        resampled=TRUE;
  } else{
        resampled=FALSE;
		}

 # Propose new particles
 
 # Chose proposal distribtion
 if(adapt & runif(1)<0.9)
 { # adapts 90% of the time
   # calculate covariance for adaptive transition kernel
   # ensure is positive definite
   proposal_cov = as.matrix(nearPD(cov(do.call(rbind,new_particles))* 2)$mat)
   # indexes correspond to the parameters we are estimating
 }
 
  proposed_particles <- lapply(particles,function(x){rmnorm(1,mean=x,varcov=proposal_cov)})      
  proposed_likelihood <- unlist(mclapply(proposed_particles,function(x)
  {
   # Test proposed values are in range 
   # two conditions- positive values and within prior ranges 
    if(min(x) < 0 | !is.finite(sum(prior_prob(x))))  
    {return(-Inf)}else
    {NBLikelihood(x,B,mu,theta,Cm,Lmax,ageGroupBreaks,StratifiedCases)}
  },mc.cores=4))

  accept_reject <- sapply(1:noParticles,function(x){
  
  # (log) Likelihood + prior + correction for log parameter transform
  current = particle_likelihood[x] + sum(prior_prob(particles[[x]])) + sum(particles[[x]])
  proposed = proposed_likelihood[x] + sum(prior_prob(proposed_particles[[x]])) + sum(proposed_particles[[x]])

   #make sure not to accept NaN values
  if(is.nan(proposed))
  {proposed=-Inf}

  if(!is.finite(proposed))  #dealing with Inf
  {compare=-Inf} else{
   # Condition should never happen
  if(current==-Inf)
  {compare=1}else
  {
    compare=exp(current-proposed);   
  }
  }

  acc=min(1,(compare));      #define the acceptance probability, scaled as appropriate
  u=runif(1);
  if(u<acc)
   {return(c(TRUE,current,proposed))}else
   {return(c(FALSE,current,proposed))} 
  })
  
  
  # Only keep accepted particles, replace rejects with particle from previous round
  proposed_particles[which(!accept_reject[1,])] = particles[which(!accept_reject[1,])]
  proposed_likelihood[which(!accept_reject[1,])] = particle_likelihood[which(!accept_reject[1,])]

 print(paste('Acceptance rate: ',mean(accept_reject)))
  
  #TransitionKernel(particleIndex)=mvnpdf(nP([2,3,4,5,9,10:12,14:15,16,17,19:20]),CurrentResampledPARTICLE([2,3,4,5,9,10:12,14:15,16,17,19:20]),TransitionKernelCovariance);
   
   TransitionKernel = sapply(1:noParticles,function(x){dmnorm(proposed_particles[[x]],mean=particles[[x]],varcov=proposal_cov)})   
  
  # Calculate new particle weights 
    if(resample)
    {new_weights=(accept_reject[2,]+abs(min(accept_reject[2,]))+1)/sum(TransitionKernel)}else{
     new_weights= weights*((accept_reject[2,]+abs(min(accept_reject[2,]))+1)/sum(TransitionKernel))}

return(list(proposed_particles,proposed_likelihood,new_weights))
     
}

