% --- Persistent Na+
% Magistretti J, Alonso A (1999) Biophysical properties and slow voltage-dependent
%   inactivation of a sustained sodium current in entorhinal cortex layer-II principal neurons. 
%   A whole-cell and single-channel study. J Gen Physiol 114:491–509
% Comment : shifted -10 mv to correct for junction potential,
%       corrected rates using q10 = 2.3, target temperature 34, orginal 21
% The voltage must be enter point by point not as a vector in order to
% eliminate singularities

function [moo,taum, hoo,tauh]=NapKinetics(V)

alpham = 0.182.*(V + 38) ./ (1 - exp(-(V + 38)./6));
betam = -0.124 .* (V + 38) ./ (1 - exp((V + 38)./6));

if abs(38 + V) < eps
  alpham = 0.182 .* 6 ./ exp(-(V + 38)./6);
  betam = 0.124 .* 6 ./ exp((V + 38)./6);
end
Tadj=(2.3).^((34-21)./10);

%moo=alpham ./ (alpham + betam);
moo=1./(1+exp(-(V + 52.6)./4.6));
%taum=1e-3 ./ (alpham + betam); % The 1e-3 factor is to transform from ms 2 s
taum = 1e-3.*6./ (Tadj .* (alpham + betam));

alphah = -2.88e-6 .* (V + 17) ./ (1 - exp((V + 17) ./ 4.63));
betah = 6.94e-6 .* (V + 64.4) ./ (1 - exp(- (V + 64.4) ./ 2.63));

if abs(17 + V) < eps
 alphah = 2.88e-6 .* 4.63 ./ exp((V + 17) ./ 4.63);
end

if abs(64.4 + V) < eps
 betah = 6.94e-6 .* 2.63 ./ exp(- (V + 64.4) ./ 2.63);
end

%hoo=alphah ./ (alphah + betah);
hoo=1./(1+exp((V + 48.8)./10));
%tauh=1e-3 ./ (alphah + betah); % The 1e-3 factor is to transform from ms 2 s
tauh =1e-3 ./ (Tadj .* (alphah + betah));

end
% to test the kinetics
 % V=-100:0.1:100; for i=1:length(V) [moo(i),taum(i), hoo(i),tauh(i)]=NapKinetics(V(i)); end; figure(1); plot(V,moo,'r.'); hold on; plot(V,hoo,'b.'); figure (2); plot(V,taum,'r.'); hold on; plot(V,tauh,'b.');
 