% Potassium delayed-rectifier -> Hodgkin and Huxley 1952
% The voltage must be enter point by point not as a vector in order to
% eliminate singularities

function [moo,taum]=KKinetics(V)
    
if V == -55
    alpham = 0.1;
else
    alpham = 0.01*(V+55)./((1-exp(-(V+55)./10)));
end

betam = 0.125*exp(-(V+65)./80);

moo=alpham ./ (alpham + betam);
taum=1e-3 ./ (alpham + betam); % The 1e-3 factor is to transform from ms 2 s

end
% to test the kinetics
% V=-180:0.1:100; for i=1:length(V) [moo(i),taum(i)]=KKinetics(V(i)); end; figure(1); plot(V,moo,'r'); figure (2); plot(V,taum,'r');
