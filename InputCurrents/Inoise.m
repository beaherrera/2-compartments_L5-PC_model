function I = Inoise(dt, mu, sigma, tau, Ipre)
% Inoise -> generates a noisy current that is Gaussian distributed with  mean mu, standard deviation sigma, and a correlation length tau 
% Inputs:
%   dt: time step
%   mu: mean of the noisy current
%   muStart: initial value of mu
%   sigma: standard deviation
%   tau: correlation length 
%   Ipre: current at t-dt, previos time


I = Ipre + (mu - Ipre).*dt./tau + sigma.*randn.*sqrt(2.*dt./tau); 

end

% to test the function
% Ipre=0; for ii=1:100000 I(ii)=Inoise(0.1, 400, 400, 3, Ipre); Ipre = I(ii); end; plot(I)