function [ mu ] = MakeMu(Lmax  )
%MakeMU generates mu for German population using Gompertz model
%INPUTS
%Lmax maximum age index

%OUTPUTS
%mu= deathr ate per age



%upper age=120
a=1:120;
DE(1)=3.232e-05;        %for 2015
DE(2)=0.09419;
mu_whole=DE(1)*exp(DE(2)*a);
mu=mu_whole(1:Lmax-1)/365;
mu(Lmax)=(-log(mean(exp(-cumsum(mu_whole(Lmax:end)))))-sum(mu_whole(1:Lmax-1))) /365;

end

