% Lytton and Sejnowski. Simulations of Cortical Pyramidal Neurons Synchronized by Inhibitory Interneurons. 1991
% Calcium dynamic is based on Kay and Wong 1987 (Hippocampal PCs)

function [moo,taum]=CaLKineticsW(V)
% Kinetics of the Ca L-type channel

alpham = 1.6./(exp(-0.072.*(V - 5))+1);
betam = 0.02.*(V + 8.69)./(exp((V + 8.69)./5.36)-1);
if abs(V + 8.69) < eps
  betam = 0.02.*5.36 ./ exp((V + 8.69)./5.36);  
end

moo=alpham ./ (alpham + betam);
taum=1e-3 ./ (alpham + betam); % The 1e-3 factor is to transform from ms->s

end
% to test the kinetics
% V=-100:0.1:130; [moo,taum]=CaLKineticsW(V); figure(1); hold on; plot(V,moo,'r'); figure(2); hold on; plot(V,taum,'r'); ECa=135; ICa=(moo.^2).*(V-ECa); figure(3); hold on; plot(V,ICa,'b');% %% Calcium L current
