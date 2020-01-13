% ---- Muscarinic K+ current, Im
%  corrected rates using q10 = 2.3, target temperature 34, orginal 21
% Adams PR, Brown DA, Constanti A (1982) M-Currents and Other Potassium 
%   Currents in Bullfrog Sympathetic Neurons. J Physiol 330: 537–572.

function [moo,taum]=ImKinetics(V)

T_adj = 2.3^((34-21)/10);

alpham = 0.0033.*exp(0.1 .* (V + 35));
betam = 0.0033.* exp(-0.1 .* (V + 35));

moo=alpham ./ (alpham + betam);
taum=1e-3./ (T_adj.*(alpham + betam));

end
% to test the kinetics
% V=-80:0.1:0; [moo,taum]=ImKinetics(V); figure(1); plot(V,moo,'r'); figure (2); plot(V,taum,'r');