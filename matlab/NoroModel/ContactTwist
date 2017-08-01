function [ TwistedContactMatrix ] = ContactTwist( ContactAgesPerAgeGroup,T,k,d)
%%
%Matrix to reduce or enhance mixing between subgroups off diagonal
%direction is as per Meyer and Held 2016
%ContacTwist takes raw contact matrix and emphasises or lessens the effect
%of within age group mixing. It also uses a mixture model for symmetric and
%assymetric contacts

%Inputs:
%ContactAgesPerAgeGroup=raw contact matrix output from MakeContactMatrix
%T=number of participants in each age group output from MakeContactMatrix
%k=Degree of diagonalisation of contact matrix, k=0 suggest C=I, k=1 suggests C is unchanged
%d=Degree of assymetry. If d=1, matrix is assymmetic, else d=0 suggests matrix is symmetric

%Output:
%TwistedContactMatrix= Partially symmetrised, partially diagonalised
%                      contact matrix
%%
noAgeGroups=length(T);

rowSums=sum(ContactAgesPerAgeGroup,2);
NORMcontactMatrix=bsxfun(@rdivide,ContactAgesPerAgeGroup,rowSums);  %normalise rows

[EiVec,EiVal]=eig(NORMcontactMatrix);   %find eigenvalues and eigenvectors of contact matrix

NewContactMatrix= (EiVec * EiVal^k) / EiVec ;  %Twist or diagonalise matrix

%truncate at zero to ensure no negative entries
NewContactMatrix(NewContactMatrix<0)=0;

%make sure matrix is real
NewContactMatrix=real(NewContactMatrix);

%un normalize or multiply by relative number of contacts
NewContactMatrix=bsxfun(@times,NewContactMatrix,rowSums);



%reciprocity- either assume no symmetry or reciprocity or complete symmetry
%or reciprocity
%DECLARE both matrix cases
AsymContactMatrix=zeros(noAgeGroups);
SymContactMatrix=zeros(noAgeGroups);
for index = 1 : noAgeGroups
    for jndex = 1 : noAgeGroups
        
        AsymContactMatrix(index,jndex) =NewContactMatrix(index,jndex)/ T(index);
        SymContactMatrix(index,jndex) =0.5 * ( NewContactMatrix(index,jndex)/ T(index) + NewContactMatrix(jndex,index)/ T(jndex) );
    end
end

%mixture model for degree of assymmetry in contacts
TwistedContactMatrix= d*AsymContactMatrix + (1-d)* SymContactMatrix;
end

