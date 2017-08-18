function [ ESS ] = ESS( weight )
%%
%ESS calculates effectivesample size based on weights

%INPUTS
%weight=vector of particle weightsbefore resampling

%OUTPUTS
%ESS=effective samplesize

%%
weight = weight / sum(weight);

ESS=1/sum(weight.^2);
end

