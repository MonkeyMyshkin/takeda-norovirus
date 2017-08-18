function [ LL ] = NBLikelihood( StratifiedCases, ReportedInfectionNumber, Dispersion)
%%
%NBLikelihood Calculates the log negative binomial likelihood for German Case notifications

%INPUTS
%StratifiedCases= German case notifications divided into same agegroups as
%                 model output
%ReportedInfectionNumber=Model output with reporting model applied
%Dispersion=Negative binomial dispersion parameter

%OUTPUTS
%LL=Log likelihood

%%
%Negative binomial function: the model cases must be in integers
%reformulate for Matlab function
R=1/Dispersion;

P=bsxfun(@rdivide,bsxfun(@times,Dispersion,StratifiedCases),bsxfun(@plus,1,bsxfun(@times,Dispersion,StratifiedCases)) );

LogLike=log( nbinpdf( round( ReportedInfectionNumber ) , R, 1-P ));

LL=sum(LogLike(:));


end


