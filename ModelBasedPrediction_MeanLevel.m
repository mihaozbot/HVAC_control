function [ymH]= ModelBasedPrediction_MeanLevel(Phi,Th,H)
%
%   function [ymH]= ModelBasedPrediction(Phi,Th,H)
%
%   Phi = [-y_z u]     regressor
%   Th                model parameters 
%   H                 število korakov predikcije
%
%   ymH               predikcija izhoda modela za H korakov



for i=1:H
    ym = Phi * Th;
    Phi = [ym Phi(2:end)];
end
ymH=ym;
