function [lfp, lfpT] = lfp_calculation(tf, Ih, ntrials, fileL5PCS)
%% Computation of the local field potentials based on the point-source approximation
% Input:
%       tf: [ms]simulation time.
%       Ih: indicates if you are processing the data generated with Ih
%           present if equal to 1, or with Ih blocked if equal to 0. This
%           variable is used to save the data.
%       ntrials: number of trials
%       Ysd: [channels x time pts] local field potential resulted from the
%           average over trials (10 trials in the case of our paper)
%       fileL5PCS: name of the folder where the data will be saved. The
%           folder will be inside the folder GenerateFigures\FiguresData
% Output:
%       lfp: [channels x time pts] local field potential resulted from the
%           average over trials (10 trials in the case of our paper)
%       lfpT: [channels X time pts X trials] local field potential of each
%           trial.
% Authors: Herrera et al., 2020

%% running trials

for trial=1:ntrials
    % Simulation time
    t0 = 0; %[ms] start time
    Fs = 10e3; % [Hz] sampling frequency
    dt = (1/Fs)*1e3;  %[ms] time increment
    
    tspan = t0:dt:(tf*1e3);   %[ms] time span
    tspan = tspan(1:(numel(tspan)-1));
    
    %% loading the simulated data
    
    cd ../
    file = ['GenerateFigures\FiguresData\' fileL5PCS];
    
    if Ih
        load([file '\currentsTrial' num2str(trial) '.mat'], 'Ysd');
    else
        load([file '\currentsNoIhTrial' num2str(trial) '.mat'], 'Ysd');
    end
    
    cd L5PCs
    
    %% compartments location
    
    Ids = Ysd; % transmembrane currents of each point-source per neuron
    ls = 890; % [um]length of the apical dendrites
    sources_number = length(Ids(1,:)); % number of sources in the column (all neurons)
    n = 5; % number of point sources per neuron
    unit_number = sources_number/n; % number of neurons
    rng default; % For reproducibility
    zc = zeros(1, sources_number); % [m] currents position, z-direction
    
    % layer 5
    % - z (um)
    a = 1025;
    b = 1450;
    rNum = rand(1,unit_number); % generating the random numbers for the z-direction
    ii = 1;
    AISZ = (b-a).*rNum + a; % this calculates the z-coordinate of the axon initial segment of each neuron
    zc(1,ii:unit_number) = (1000-700).*rNum + 700; % this calculates the z-coordinate of the soma
    zc(1,ii*unit_number+1:(ii+1)*unit_number) = AISZ; ii = ii+1; % initial segment source position
    zc(1,ii*unit_number+1:(ii+1)*unit_number) = AISZ + 150; ii = ii+1; % basal compartment source position, 150 um below the axon initial segment
    zc(1,ii*unit_number+1:(ii+1)*unit_number) = AISZ - ls; ii = ii+1; % main-bifurcation point, ls um above the soma
    zc(1,ii*unit_number+1:(ii+1)*unit_number) = AISZ - ls - 150; % tuft, 150 um above the main-bifurcation point
    
    % -- calculating the position of each neuron on the plane. The same
    % coordinates are used for all sources per neuron
    % - x and y coordinates (mm)
    rc = 3/2; % [mm] radius of the cortical column
    xc = zeros(1,unit_number);
    yc = zeros(1,unit_number);
    
    for ii=1:unit_number
         r = sqrt(rand())*rc;
         theta = rand()*2*pi;
         xc(ii) = r*sin(theta);
         yc(ii) = r*cos(theta);
    end
    
    xc = repmat(xc,1,n);
    yc = repmat(yc,1,n);
    
    % converting to m
    zc = zc.*1e-6; % um-> m
    yc = yc.*1e-3; % mm -> m
    xc = xc.*1e-3; % mm -> m
    
    
    %% Electrodes position
    Ne = 16; % number of electrodes in the shank
    a = 0.1; % [mm] position of the first electrode
    elec_spacing = 0.1; % [mm] electrode spacing
    ze = a:elec_spacing:((Ne-1)*elec_spacing + a); % position of the electrodes along z with respect to pia matter
    el_pos = ze*1e-3;  % [m] electrode positions with respect to the pia surface
    cond = 0.323; %[S/m] gray matter conductance
    
    %% Current sources
    Iz = Ids';
    
    %% LFP calculation
    if trial==1
        lfpT = zeros(length(ze), length(Iz(1,:)), ntrials); % preallocating the local field potential matrix
    end
    
    lfpt = zeros(length(ze), length(Iz(1,:)));
    % Emapping = zeros(length(ze), sources_number);
    
    for ii=1:Ne
        [lfpt(ii,:), ~] = lfp_pointSourceMethod(Iz, xc, yc, zc, cond, el_pos(ii), elec_spacing.*1e-3); % mV
    end
    
    lfpT(:,:,trial) = lfpt;
    
end

%% averaging the LFPs across trials
lfp = mean(lfpT, 3);

%% visualization

font = 16;

figure;
Ve = eegfilt(lfp, Fs, 0, 90,0,250); % filters the LFP at 90 Hz

plot(tspan, Ve.*1e3);
xlabel('Time (ms)', 'FontSize',font);
ylabel('LFP (\muV)','FontSize',font);
title(['Mean LFP across ' num2str(ntrials) ' trials'])

%% save

cd ../
file = ['GenerateFigures\FiguresData\' fileL5PCS];
if ~exist(file, 'dir') % checks if the folder already exists
    mkdir(file);  % creates a folder named 'file'
end

if Ih % Ih present
    save([file '\lfp_L5PCs.mat'], 'lfp');
    save([file '\lfp_L5PCsTrials.mat'], 'lfpT');
else % Ih blocked
    save([file '\lfp_L5PCsNoIh.mat'], 'lfp');
    save([file '\lfp_L5PCsNoIhTrials.mat'], 'lfpT');
end

cd L5PCs

end
