require(rootSolve)
require(deSolve)
require(Rcpp)

sourceCpp('../cpp/NoroModel.cc')

AgeStratify <- function(Cases, ageGroupBreaks)
{

# AgeStratify groups cases into by age group from annual counts
# Inputs: Cases = a vector or matrix of case numbers for each annual age
#                group
#        ageGroupBreaks= a vector of age group divisions 
# Outputs: StratifiedCases = counts grouped by age group
# Translated from matlab code (almost certainly a more elegant solution)

# Trap for vector input
if(is.null(dim(Cases)))
{ 

 StratifiedCases = matrix(length(ageGroupBreaks))

 maxAge=length(Cases);

# need upper bound on ages as limit for collecting counts
ageGroupBreaks = c(ageGroupBreaks, maxAge)

# Calculate the number of age groups
noAgeGroups=length(ageGroupBreaks)

#for each week, group cases into each age group ignoring NaN
StratifiedCases[1]= sum(Cases[1 : ageGroupBreaks[1]])

for(index in c(2 : (noAgeGroups-1)))
{
    StratifiedCases[index]= sum(Cases[(ageGroupBreaks[index-1]+1) : (ageGroupBreaks[index])]);

}

return(StratifiedCases)

} else{ 

StratifiedCases = matrix(NA,dim(Cases)[1],length(ageGroupBreaks))

maxAge=dim(Cases)[2]

# need upper bound on ages as limit for collecting counts
ageGroupBreaks = c(ageGroupBreaks, maxAge)

# Calculate the number of age groups
noAgeGroups=length(ageGroupBreaks)

#for each week, group cases into each age group ignoring NaN
StratifiedCases[,1]= rowSums(Cases[, 1 : ageGroupBreaks[1]])


for(index in c(2 : (noAgeGroups-1)))
{
    StratifiedCases[,index]= rowSums(Cases[, (ageGroupBreaks[index-1]+1) : (ageGroupBreaks[index])]);

}

return(StratifiedCases)

}


}

MakeContactMatrix <- function(ageGroupBreaks,GermanContactData,GermanParticipantData)
{
#
# MakeContactMatrix takes Polymod data and produces contact matrix and age
# group counts
#
#Inputs:
#(POLYMOD data: global ids from contacts and participants, participant and
#contact ages and weekend/weekday info: weekday=1, weekend=0)
#GermanContactData=[GlobalIDs, ContactAges]
#GermanParticipantData=[GlobalIds, ParticipantAges , Weekday?]
#Outputs:
#
# returns list with three elements:
#contactMatrix= a symmetrised average probability of contact between age
#               groups see Bageulin et al. 2013
# ContactAgesPerAgeGroup=a raw matrix of contacts between each age group
#T= Total numbers of participants in each age group

#define a maximum age for age counts
maxAge=99;

ageGroupBreaks=c(ageGroupBreaks, maxAge)

noAgeGroups=length(ageGroupBreaks)

noParticipants=length(GermanParticipantData[,1])

# sort participant ages

GermanParticipantData= GermanParticipantData[sort.int(GermanParticipantData[,2],index=T)$ix,]


# for each participant, we note whether they have contact with a contact of
# each yearly age
# DECLARE: row for every participant, column for every age

ContactAgesPerParticipant=matrix(0,noParticipants,noAgeGroups);

for(index in c(1 : noParticipants))
{
    #identify all contacts related to current participant
    #isolate the ages of these contacts
    #sort out frequency of contacts in each age group: agegroupbreaks are
    #increased by 1 as they correspond to an index initialised from 1
    ContactAgesPerParticipant[index,] = hist(GermanContactData[GermanContactData[,1]==GermanParticipantData[index,1],2],c(-1,ageGroupBreaks),plot=FALSE)$counts;
    

}

#Define the number of participants of each age, T, in each age
#group, W, and number of weekday surveys, Nwd

W=hist(GermanParticipantData[,2],c(0:(maxAge+1))-0.5,plot=F)$counts
T=hist(GermanParticipantData[,2],c(0,ageGroupBreaks+1)-0.5,plot=F)$counts
Nwd=sum(GermanParticipantData[,3]!=0)   #number of non zero elements


# Calculate weightings for weekday/weekend participant surveys
# DECLARE: weight diagonal matrix for ease later on

weight=matrix(0,noParticipants,noParticipants)

for (index in c(1:noParticipants))
{
    
    ageGroupIndex= which(ageGroupBreaks>GermanParticipantData[index,2])[1];
    
    if(GermanParticipantData[index,3])
    {
      weight[index,index] = (W[GermanParticipantData[index,2]+1]/T[ageGroupIndex] ) * 5/Nwd
    }else{
        weight[index,index] = (W[GermanParticipantData[index,2]+1]/T[ageGroupIndex] )* 2/(noParticipants-Nwd)
    }
    
}

#Update ContactAgesPerParticipant with weighting (matrix multiplication)

ContactAgesPerParticipant=weight %*% ContactAgesPerParticipant;


# Calculate all contacts for participants in one age group with contacts in
# another
# DECLARE: row , column for every age

ContactAgesPerAgeGroup=matrix(0,noAgeGroups,noAgeGroups)

# first need indexes of last participants in each age group

partIndexLast=numeric(noAgeGroups)

for(index in c(1 : noAgeGroups))
{
    partIndexLast[index]=rev(which(GermanParticipantData[,2]<=ageGroupBreaks[index]))[1];
}


ContactAgesPerAgeGroup[1,] = colSums( ContactAgesPerParticipant[ c(1 : partIndexLast[1]), ]) / sum(diag( weight[ 1 : partIndexLast[1], c(1 : partIndexLast[1])]));


for(index in c(2:noAgeGroups))
{

   ContactAgesPerAgeGroup[index,] = colSums( ContactAgesPerParticipant[ c((partIndexLast[index-1]+1) : partIndexLast[index]),]) / sum(diag( weight[c(( partIndexLast[index-1]+1):partIndexLast[index]),c((partIndexLast[index-1]+1): partIndexLast[index])]))
    
}

#finally symmetrise contactmatrix
#Declare: row, column for every age

contactMatrix=matrix(0,noAgeGroups,noAgeGroups);

for(index in c(1:noAgeGroups))
{
    for(jndex in c(1:noAgeGroups))
    {
        
        contactMatrix[index,jndex] = ( ContactAgesPerAgeGroup[index,jndex] * T[index] + ContactAgesPerAgeGroup[jndex,index]* T[jndex] )/(2*T[index]);
        
    }
}

return(list(contactMatrix,ContactAgesPerAgeGroup,T))

}

MakeMu <- function(Lmax)
{
#MakeMU generates mu for German population using Gompertz model
#INPUTS
#Lmax maximum age index

#OUTPUTS
#mu= death rate per age
#upper age=120
a=1:120
DE = c(3.232e-05,0.09419)  #for 2015

mu_whole=DE[1]*exp(DE[2]*a)
mu = mu_whole[1:Lmax-1]/365
mu[Lmax]=(-log(mean(exp(-cumsum(mu_whole[Lmax:120]))))-sum(mu_whole[1:(Lmax-1)]))/365

return(mu)

}

ContactTwist <- function(ContactAgesPerAgeGroup,T,k,d)
{

# Matrix to reduce or enhance mixing between subgroups off diagonal
# direction is as per Meyer and Held 2016
# ContacTwist takes raw contact matrix and emphasises or lessens the effect
# of within age group mixing. It also uses a mixture model for symmetric and
# asymmetric contacts

# Inputs:
# ContactAgesPerAgeGroup=raw contact matrix output from MakeContactMatrix
# T=number of participants in each age group output from MakeContactMatrix
# k=Degree of diagonalisation of contact matrix, k=0 suggest C=I, k=1 suggests C is 
# unchanged
# d=Degree of assymetry. If d=1, matrix is asymmetric, else d=0 suggests matrix is 
# symmetric

# Output:
# TwistedContactMatrix = Partially symmetrised, partially diagonalised
#                       contact matrix
#

# TwistedContactMatrix

noAgeGroups=length(T)

NORMcontactMatrix= ContactAgesPerAgeGroup / rowSums(ContactAgesPerAgeGroup)

#find eigenvalues and eigenvectors of contact matrix
oot <- eigen(NORMcontactMatrix)

EiVec <- oot$vectors
EiVal <- diag(oot$values)



NewContactMatrix= (EiVec * EiVal^k) %*% solve(EiVec)  #Twist or diagonalise matrix


# truncate at zero to ensure no negative entries
NewContactMatrix[NewContactMatrix<0]=0;

#make sure matrix is real
NewContactMatrix=Re(NewContactMatrix);

#un normalize or multiply by relative number of contacts
NewContactMatrix=NewContactMatrix * rowSums(ContactAgesPerAgeGroup)




# reciprocity- either assume no symmetry or reciprocity or complete symmetry
#or reciprocity
# DECLARE both matrix cases
AsymContactMatrix=matrix(0,noAgeGroups,noAgeGroups)
SymContactMatrix=matrix(0,noAgeGroups,noAgeGroups)

for(index in c(1 : noAgeGroups))
{
    for(jndex in c(1 : noAgeGroups))
    {   
        AsymContactMatrix[index,jndex] = NewContactMatrix[index,jndex];
        SymContactMatrix[index,jndex] =( ContactAgesPerAgeGroup[index,jndex]* T[index] + ContactAgesPerAgeGroup[jndex,index]* T[jndex] )/(2*T[index]);
    }
}

#mixture model for degree of assymmetry in contacts
TwistedContactMatrix= d*AsymContactMatrix + (1-d)* SymContactMatrix;
return(TwistedContactMatrix)
}

# Helper functions to pack 2d matrix into array
index2 <- function(i,a,y)
{
  return (y*(i-1) + a)
}



Expand <- function(Cm,ageGroupBreaks,max_age)
{

# Function to expand coarse-grained age matrix Cm
# to annual age-cohorts

# Input: Cm, contact matrix with age groups defined by
#        ageGroupBreaks (vector)
#        max_age (numeric)

# Output: Matrix with annual age groups ( max_age, max_age)

CmEx <- matrix(0,max_age,max_age)

ageGroupBreaks <- c(1,ageGroupBreaks,max_age)
	
	for(i in 1:(length(ageGroupBreaks)-1) )
	{
	 for(j in 1:(length(ageGroupBreaks)-1))
	 {
	  for(a1 in seq(ageGroupBreaks[i],ageGroupBreaks[i+1],1))
	  {
	  for(a2 in seq(ageGroupBreaks[j],ageGroupBreaks[j+1],1))
	    {
	     CmEx[a1,a2] = Cm[i,j];
	    }
	  }
     }
	}
	return(CmEx)
}
	

mk_initial_conditions <- function(max_age,params,Cm,mu,B,theta)
{

a = 1/365.0

x = numeric(7*max_age)

for(k in c(1:7))
{
x[index2(k,1,max_age)] = (B/(mu[1] + a))/7.0;
}
  
for( a1 in c(2:(max_age-1)))
{
     for(k in c(1:7))
    	{
    	 x[index2(k,a1,max_age)] = ((a*x[index2(k,a1-1,max_age)])/(mu[a1]+a));
     	}
}
 
for(k in c(1:7))
{
	x[index2(k,max_age,max_age)] = (a*x[index2(k,max_age-1,max_age)]/(mu[max_age]));
}

FPs<-steady(x,0,NoroFP,parms=params,positive=TRUE,Cm=Cm,mu=mu,max_age=max_age,B=B,theta=theta)$y

dim(FPs) <- c(max_age,7)

return(cbind(FPs,0,rowSums(FPs)))
}

              

SimulateSeasons <- function(params, B, mu,theta,ContactMatrix,max_age,burnin)
{
# SimulateSeasons: Takes input parameters and initial conditions and solves the system of 
# ODEs for all
# seasons and all age groups

# INPUTS:
# Param=vector of parameter values. Doubles.
# mu=vector of death rates for each age group. Doubles.
# theta=vector of susceptibility of recovered individuals for each season.
#       Doubles.
# ContactMatrix= noAgeGroups x noAgeGroups matrix of contact probabilities
#               between age groups
# x0=vector of initial conditions for each disease state
# omega2= vector of seasonal offsets for each age group

# OUTPUTS:
# SimulationResult=solution of ODEs for each week in TIME and each age
#

init <- mk_initial_conditions(max_age,params,ContactMatrix,mu,B,0.0)
x <- as.numeric(init)

# Forward simulate for burnin*52 weeks + 19 to bring up to start of reported data
oot<-try(lsoda(x,seq(0,burnin*365,7.0),NoroDyn,parms=params,Cm=ContactMatrix,mu=mu,max_age=max_age,B=B,theta=0))

# Season 1

oot<-try(lsoda(oot[dim(oot)[1],2:dim(oot)[2]],seq(0,33*7,7),NoroDyn,parms=params,Cm=ContactMatrix,mu=mu,max_age=max_age,B=B,theta=theta[1]))

cases = diff(oot[,1+index2(9,1:81,max_age)])

# Seasons 2:8

for(i in 2:8)
{

oot<-try(lsoda(oot[dim(oot)[1],2:dim(oot)[2]],seq(33*7+(52*7)*(i-2),33*7 +(52*7)*(i-1),7),NoroDyn,parms=params,Cm=ContactMatrix,mu=mu,max_age=max_age,B=B,theta=theta[i]))

cases = rbind(cases,diff(oot[,1+index2(9,1:81,max_age)]))

}

# Season 9

oot<-try(lsoda(oot[dim(oot)[1],2:dim(oot)[2]],seq(33*7 +(52*7)*(8-1),33*7 +(52*7)*(8-1) + 21*7,7),NoroDyn,parms=params,Cm=ContactMatrix,mu=mu,max_age=max_age,B=B,theta=theta[9]))

cases = rbind(cases,diff(oot[,1+index2(9,1:81,max_age)]))

#if(is.null(dim(oot)) || dim(oot)[1] != 209 ){return(NA)}       
#else{ cases = oot[,seq(1,81)*7 + 1]}        

return(cases)

}

