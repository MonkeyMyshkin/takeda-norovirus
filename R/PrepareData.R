source('ModelFunctions.R')

# Data

GermanCaseNotification<-read.csv('../data/GermanCaseNotification.csv',header=FALSE)
GermanPopulation <- scan('../data/GermanPopulation.csv')
GermanContactData<-read.csv('../data/GermanContactData.csv',header=FALSE)
GermanParticipantData<-read.csv('../data/GermanParticipantData.csv',header=FALSE)

#AGE GROUPS
ageGroupBreaks=c(8,18,26,37,50,70)

noAgeGroups=length(ageGroupBreaks)+1
#calculate population sizes for these age groups
PopulationSize=AgeStratify(GermanPopulation,ageGroupBreaks)
StratifiedCases<-AgeStratify( GermanCaseNotification, ageGroupBreaks )

#CONTACT MATRIX

oot <- MakeContactMatrix( ageGroupBreaks, GermanContactData , GermanParticipantData )

contactMatrix = oot[[1]]
ContactAgesPerAgeGroup = oot[[2]]
T = oot[[3]]
rm(oot)

#DEATH FUNCTION
#maximum age group in case notifications is 80. Therefore our max age will
#be 80, indexed as 81.
MaxAge=80;
Lmax=MaxAge+1;
mu=MakeMu(Lmax);    #daily death rate

# Expand contact matrix to model age-structure
Cm = Expand(contactMatrix,ageGroupBreaks,Lmax)
       
B = (17.0/1000)/365.0

    
        
