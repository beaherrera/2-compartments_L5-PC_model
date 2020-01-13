% --- Non-specific cation current
% Kole MH, Hallermann S, Stuart GJ (2006) Single Ih channels in pyramidal 
%   neuron dendrites: properties, distribution, and impact on action potential
%   output. J Neurosci 26: 1677–1687.
% The voltage must be enter point by point not as a vector in order to
% eliminate singularities

function[moo,taum]=IhKinetics(V)

alpham = (6.43 .* (V + 154))./ (exp((V + 154)./11.9)-1);
betam = 193.*exp(V./33.1);

if isnan(alpham)
    alpham = (6.43)./ ((1./11.9).*exp((V + 154)./11.9));
end

moo= alpham./(alpham + betam);
taum=1./(alpham + betam);

end
% to test the kinetics
 % V=-220:0.1:120; for i=1:length(V) [moo(i),taum(i)]=IhKinetics(V(i)); end; figure(1); plot(V,moo,'r'); figure (2); plot(V,taum,'r');
 