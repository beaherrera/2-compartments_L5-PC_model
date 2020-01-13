% --- Slow inactivating K+ current
% Comment : shifted -10 mv to correct for junction potential,
%       corrected rates using q10 = 2.3, target temperature 34, orginal 21.
% Korngreen A, Sakmann B (2000) Voltage-gated K+ channels in layer 5 
%    neocortical pyramidal neurones from young rats: subtypes and gradients. 
%    J Physiol 525: 621–639.
% The voltage must be enter point by point not as a vector in order to
% eliminate singularities

function [moo,taum, hoo, tauh]=KslowKinetics(V)

Tadj=(2.3).^((34-21)./10);

moo=1./(1+exp(-(V + 11)./12));
hoo=1./(1+exp((V + 64)./11));

if V < -50
    taum=1e-3 .* (1.25 + 175.03 .* exp(0.026 .* (V + 10))) ./ Tadj; % The 1e-3 factor is to transform from ms 2 s
else    
    taum=1e-3 .* (1.25 + 13 .* exp(-0.026 .* (V + 10))) ./ Tadj;   
end

tauh=1e-3 .* (360 + (1010 + 24 .* (V + 65)).* exp(-((V + 85)./48).^2)) ./ Tadj;

end
% to test the kinetics
% V=-100:0.1:100; for i=1:length(V), [moo(i),taum(i), hoo(i),tauh(i)]=KslowKinetics(V(i)+20); end; figure(1); plot(V,moo,'r.'); hold on; plot(V,hoo,'b.'); figure (2); plot(V,taum,'r.'); hold on; plot(V,tauh,'b.');
