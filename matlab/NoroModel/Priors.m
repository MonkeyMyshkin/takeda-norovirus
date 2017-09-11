function [ PriorProbabilities ] = Priors(Param,theta,ReportingParam,Dispersion,k,d,damping,GermanPopulation,ContactMatrix,mu,ageGroupBreaks)
%%PRIORS takes proposed values of each parameter and returns sum of log prior probability

%INPUTS
%Param=vector including alpha, q,omega1,omega2,nu, delta, epsilon, sigma
%               and psi 
%theta=vector of susceptibility of recovered individuals per season
%ReportingParam= reporting baseline values with one value for 18-26 and
%                26-37
%Dispersion=Dispersion parameter for negative binomial likelihood
%k=degree of diagonalisation of contact matrix
%d=degree of assymmetry in contact matrix
%damping= damping(1) is scaling on damping for 0-37 year olds and
%         damping(2) is baseline reporting dampin factor- see CappedReporting
%GermanPopulation=Age group sizes for German Population
%ContactMatrix= contact matrix given current parameters
%mu=vector of death rates per age
%ageGroupBreaks=vectorof age group divisions

%OUTPUTS
%PriorProbabilities=sum of log prior probabilities 
%%
%seasonal forcing amplitude
if Param(3)>0.01 && Param(3)<0.3
    omega1Prior = log(gampdf(Param(3),0.15,1)/(gamcdf(0.3,0.15,1)-gamcdf(0.01,0.15,1)));
else
    omega1Prior = log(0);
end


%scaling of asymptomatic infectiousness
if Param(4)>0 && Param(4)<1
    nuPrior = log(gampdf(Param(4),1,1)/ (gamcdf(1,1,1)-gamcdf(0,1,1)));
else
    nuPrior=log(0);
end

%loss of immunity
if  Param(5)>1/(10*365) && Param(5)< 2/(365) %longer than 6 months, less than ten years
    deltaPrior = log( gampdf(Param(5),1/(2*365),1)/ (gamcdf(2/365,1/(2*365),1)-gamcdf(1/(10*365),1/(2*365),1)));
else
    deltaPrior =log(0);
end
%deltaPrior = log( unifpdf(ParamVector(6),1/(30*365),1));

%proportion of symptomatic
% if Param(7)<=1 && Param(7)>=0
%     sigmaPrior = log(normpdf(Param(7),0.735415,0.0960897));    %from challenge study data and ASYMP_PROB program
% else
%     sigmaPrior=log(0);
% end
% sigmaPrior = log(unifpdf(ParamVector(8)));


%recovery rate
gammaPrior = log( gampdf(Param(9),1,1));

%ReportingParam
ReportingPrior(1:6)=log(unifpdf(ReportingParam,0,1));  %uniform prior just to keep within feasible limits
    %prior on average reporting rate
ReportingPrior(7)=log(unifpdf(mean(ReportingParam(1)*[1 ReportingParam(2)*ones(1,2) ReportingParam(3) ReportingParam(4) ReportingParam(5) ReportingParam(6)]),0,0.1));


%Dispersion
DispersionPrior=log(unifpdf(Dispersion,0,0.5));


%k
kPrior=log( unifpdf(k,0,2));

%d
% dPrior=log(unifpdf(d));


%theta
% for i=1:9
%     if theta(i)<1
%         thetaPrior(i)=log(normpdf(theta(i),1,0.1)/(normcdf(1,1,0.1)-normcdf(0,1,0.1)) );
%     else
%         thetaPrior(i)=log(0);
%     end
% end

%damping
dampPrior=[log(unifpdf(damping(1),0,1)),log(gampdf(damping(2),1e-2,1e-1))];

%prior on R0
R0=MakeR0( Param(1:9) ,GermanPopulation,ContactMatrix, mu, ageGroupBreaks);
if R0>1
    R0Prior=log(normpdf(R0,15,1));
else
    R0Prior=log(0);
end

%sum of log prior probabilities
PriorProbabilities=  omega1Prior  +nuPrior+ deltaPrior + sum(ReportingPrior) ...
    +sum(dampPrior) + gammaPrior+kPrior +DispersionPrior + R0Prior;

end

