function [lfp, mapping] = lfp_pointSourceMethod(In, xc, yc, zc, sigma, el_pos, h)
%% Calculates the extracellular electric potential using the point source methods.
%% Parameters
% - Inputs:
    % In:  point current source (mA) -> (array)[point sources]x[time points]
    % xc: neuronal compartment position, x-axis (m) -> (array) [1]x[point sources]
    % yc: neuronal compartment position, y-axis (m) -> (array) [1]x[point sources]
    % zc: neuronal compartment position, z-axis (m) -> (array) [1]x[point sources]
    % sigma: extracellular conductivity (Assumption: constant through out the
    % cortical column) (S/m) (scalar)
    % el_pos: extracellular electrodes position, z-axis. (x=0,y=0) (m) (scalar)
    % h: inter-electrode distance. (m)
% - Outputs:
    % lfp: extracellular potential (mV) -> [1 channel]x[time points]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vc = 2*pi*(0.4/2)^2; % mm3

mapping = sqrt(xc.^2 + yc.^2 + (zc-el_pos).^2) - abs(zc-el_pos); 
lfp = (h/(2.*sigma)).*sum(In.*mapping')./(Vc*1e-9);

end
