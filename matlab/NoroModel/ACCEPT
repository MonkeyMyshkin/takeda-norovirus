function [ accept ] = ACCEPT(  AccProp, AccCurrent)
%%
%ACCEPT completes accept or reject step
%
%INPUTS
%AccProp=proposed likelihood and prior probabilities
%AccCurrent=current likelihood and prior probabilities

%OUTPUTS
%accept=binary term, 1 if accept and 0 if reject
%%

if isnan(AccProp)>0        %make sure not to accept NaN values
    AccProp=-Inf;
end

if AccProp==-Inf           %dealing with Inf
    compare=-Inf;
elseif AccCurrent==-Inf
    compare=1;
else
    
    compare=exp(AccProp-AccCurrent);   
end

acc=min(1,(compare));      %define the acceptance probability, scaled as appropriate

u=rand(1);

if u<acc
    accept=1;
else
    accept=0;
end
end

