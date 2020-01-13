% Na Hodgkin and Huxley 1952
% The voltage must be enter point by point not as a vector in order to
% eliminate singularities

function [moo,taum, hoo,tauh]=NaKinetics(V)
% Resting Membrane Potential = -65 mV

if V == -40
    alpham = 1;
else
    alpham = (V + 40)./(10.*(1-exp(-(V+40)./10)));
end

betam = 4*exp(-(V+65)/18);

moo=alpham ./ (alpham + betam);
taum=1e-3 ./ (alpham + betam); % The 1e-3 factor is to transform from ms 2 s

alphah = 0.07*exp(-(V+65)/20);
betah = 1./(1 + exp(-(V+35)/10));

hoo=alphah ./ (alphah + betah);
tauh=1e-3 ./ (alphah + betah); % The 1e-3 factor is to transform from ms 2 s

end

 % V=-180:0.1:120; [moo,taum, hoo,tauh]=NaKinetics(V); figure(1); plot(V,moo,'r.'); hold on; plot(V,hoo,'b.'); figure (2); plot(V,taum,'r.'); hold on; plot(V,tauh,'b.');
 