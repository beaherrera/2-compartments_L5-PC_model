%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collection of Layer 5 Pyramidal cells described by our mininal model
%
% Compartments: basal-dendritic/somatic compartment and
% apical-dendrite/trunk compartment
% - Soma/basal-dendrites compartment includes typical Na+ and K+ conductances. 
% - Trunk/apical-dendrites compartments includes: persistent Na+ (Nap), 
%   hyperpolarization-activated cation (Ih), slow inactivation K+ (Ks), 
%   muscarinic K+ (IM) and Ca2+ L-type currents (CaL). 
%
% created by Herrera et al., 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
% close all

%% Global variables
global SigX const_current pulse_train noisy_current
global I3_sq I2_sq I3pre
global dt
global I_trans ICaL
global non_state_vars non_state_num
global I2_Final I3_Final

%% simulation time
t0 = 0; %[s] start time
tf = 80*1e-3;  %[s] end time 
dt = 0.000001;  %[s] time increment
CSR=round(0.001./dt); % downsampling the tspan to export [Ca2+]i  

tspan = t0:dt:tf;   %[s] time span
tspan = tspan(1:(numel(tspan)-1));

%% stochastic or deterministic
stochastic = true;  % if true system is stochastic - if false system is deterministic

%%  ----------------------------------------------- parameters

unit_number = 10; % number of L5-PCs
ntrials = 1; % number of trials

Rd_i = 65e6; % [Ohm] transfer resistance between the Soma/basal-dendrites and Trunk/apical-dendrites compartments

I_trans = zeros(length(tspan),unit_number); % transfer current between soma and apical tuft per neuron
ICaL = zeros(length(tspan),unit_number); % Ca2+ current per neuron

rng('shuffle')

fileL5PCS = ['data' num2str(unit_number) 'n']; % name of the folder where the data will be saved

f = 149; % frequency of the train of pulses when pulse_train =1

%% -- stimulus type
const_current = 0; % if true(=1): soma stimuli is a current step pulse, and the
%                    dendritic stimuli is either a current step pulse or a double exponential
pulse_train = 0; % if true(=1): the somatic/dendritic input is a train of square current pulses
noisy_current = 1; % if true(=1): the stimuli is a noisy current with constant mean and std
% When all variables are false(=0) the somatic/dendritic input
% current is a Noisy staircase current. This type of current was used to
% generate the f-I curve for either somatic or dendritic stimulations.

%% ----------------------------------------------- initial conditions (ICs)
% ICs are set here. However after defining the parameters in the next
% section, the steady state of the system is found using fsolve and the
% results will be used as the actual ICs.

Vs_0 = -65.*ones(1, unit_number); % somatic voltage
Vd_0 = -65.*ones(1, unit_number); % dendritic voltage
[NamS_0, ~, NahS_0, ~] = NaKineticsB(Vs_0); % Na gating variables
[KmS_0,~]=KKineticsB(Vs_0); % K gating variables
[CamD_0, ~] = CaLKineticsW(Vd_0); % CaL gating variables
Cai_inf = 76.6e-6; % [mM] intracellular equilibrium calcium concentration at rest
Cai_0 = Cai_inf.*ones(1, unit_number); % [mM]
[NapmD_0, ~, NaphD_0, ~] = NapKinetics(Vd_0); % Nap gating variables
[KsmD_0, ~, KshD_0, ~] = KslowKinetics(Vd_0); % Ks gating variables
[Ihm_0, ~]=IhKinetics(Vd_0); % Ih gating variables
[Imm_0, ~]=ImKinetics(Vd_0); % Im gating variables

X0 = [Vs_0, Vd_0, NamS_0, NahS_0, KmS_0, CamD_0, Cai_0, NapmD_0, NaphD_0, KsmD_0, KshD_0, Ihm_0, Imm_0]';

state_var_num = numel(X0);  % number of state variabels

%% Variance of Wiener Process
% We considered only the intracellular calcium concentration as an
% stochastic process.

g_Vs = 0.05.*ones(1, unit_number);
g_Vd = 0.02.*ones(1, unit_number);
g_NamS = 0.*ones(1, unit_number);
g_NahS = 0.*ones(1, unit_number);
g_KmS = 0.*ones(1, unit_number);
g_CamD = 0.*ones(1, unit_number);
g_Cai = 0.000000001.*ones(1, unit_number);
g_NapmD = 0.*ones(1, unit_number);
g_NaphD = 0.*ones(1, unit_number);
g_KsmD = 0.*ones(1, unit_number);
g_KshD = 0.*ones(1, unit_number);
g_IhmD = 0.*ones(1, unit_number);
g_ImmD = 0.*ones(1, unit_number);

stdInstCa = 0.00001; % intrumental noise of the Ca2+ traces


if stochastic
    SigX = [g_Vs, g_Vd, g_NamS, g_NahS, g_KmS, g_CamD, g_Cai, g_NapmD, g_NaphD, g_KsmD, g_KshD, g_IhmD, g_ImmD]'./sqrt(dt);
else
    SigX = zeros(state_var_num,1);
end

%% running trials

for trial=1:ntrials
    
    clearvars Xs Ys Ysd
    
    sprintf('trial=%d', trial)

%% --------------------------------------------- Stimulus
% these variables store the dendritic and somatic stimulus for later
% visualization
I2_Final=zeros(size(tspan));
I3_Final=zeros(length(tspan),unit_number);

%% ------------------------ Constant Current
if const_current
    
    % -- BACfiring reproduction
    %   paradigm=1; subthreshold dendritic stimulation only
    %   paradigm=2; suprathreshold somatic stimulation only
    %   paradigm=3; combine subthreshold dendritic and suprathreshold
    %               somatic stimulations
    %   paradigm=4; suprathreshold dendritic stimulation only
    
    paradigm = 1;
    
    % --------- Soma, input current
    I3_amp = 1e-6;% [mA] amplitude of the somatic input I3; we used I3_amp=1e-6 mA in the paper
    I3_stim_onset = 30e-3; % [ms] stimulus onset
    I3_stim_end = 35e-3; % [ms] stimulus end
    % Note: the duration of the current pulse used in the paper was 5 ms.
    
    % ---------------- Apical, input current
    % The dendritic input can be a current pulse or a continuos current
    % generated by a double exponential function of the form I2_amp*(1-exp(-1/tau1))*exp(-1/tau2) 
    % with tau1=2 ms and tau2=10 ms.
    % When I2_stim_end=0, the code selects the double exponential as input
    % current.
    I2_amp = 0.6e-6; % [mA] amplitude of the dendritic input I2; threshold for Ca-spike=~0.805e-6 mA for the double exponential case
    I2_stim_onset = 38e-3; % [ms] stimulus onset
    I2_stim_end = 0; % [ms] stimulus end
    % Note: for the double exponential case the amplitude of the stimulus
    % is not the value you enter since you multiply that amplitude by the
    % double exponential function
    
    I_parameters = [I3_amp I3_stim_onset I3_stim_end I2_amp I2_stim_onset I2_stim_end]; % to enter non-state variables into the function
    
    %% ------------------------ Train, pulses
    
elseif pulse_train
    
    % We used the train of square current pulses at different frequencies
    % to study the influence of the somatic APs on the dendritic Ca2+
    % spike generation.
    % The critical frequency when Ih is blocked was 149 Hz.
    
    len = length(tspan); % total number of time points
    
    % --------- Soma, input current
    I3_amp = 15e-6; % [mA] amplitude of the somatic input current, I3
    widthAct_I3 = 2e-3; % [s] width of the pulses
    Freq_I3 = f; % [Hz] pulses frequency
    I3_stim_onset = 1e-3; % [s] stimulus onset
    I3_stim_end = 100e-3; % [s] stimulus end
    PWindowI3=round([I3_stim_onset,I3_stim_end]./dt);
    width_I3 = round(widthAct_I3./dt); % width of the square pulse in time pts
    IP_I3 = round(1./(Freq_I3.*dt) - width_I3); % inter-pulse interval
    NC_I3 = round((tf./dt)./(width_I3 + IP_I3)); % number of pulses
    
    I3_sq = square_pulse(I3_amp, width_I3, IP_I3, NC_I3, len, PWindowI3); % generates the square pulses
    
    % ---------- Apical, input current
    I2_amp = 0; % [mA] amplitude of the dendritic input current, I2
    widthAct_I2 = 2e-3; % [s] width of the pulses
    Freq_I2 = 20; % [Hz] pulses frequency
    I2_stim_onset = 30e-3; % [ms] stimulus onset
    I2_stim_end = 50e-3; % [ms] stimulus end
    PWindowI2=round([I2_stim_onset,I2_stim_end]./dt);
    width_I2 = round(widthAct_I2./dt); % width of the square pulse in time pts
    IP_I2 = round(1./(Freq_I2.*dt) - width_I2); % inter-pulse interval
    NC_I2 = round((tf./dt)./(width_I2 + IP_I2));
    
    I2_sq = square_pulse(I2_amp, width_I2, IP_I2, NC_I2, len, PWindowI2); % generates the square pulses
    
    I_parameters = [];
    
    %% ------------------------ Current with noise
    
elseif noisy_current
    
    % --------- Soma, input current
    I3_mu = 90e-6; % effective current amplitude. This amplitude is multiply by randn resulting in a noisy current with amplitudes between -5 and 5nA.
    I3_sigma = 0.2e-6; % [mA] current standard deviation
    I3_tau = 3e-3; % [s] correlation length
    I3_onset = 10e-3; % [s] stimulus onset
    I3_end = 30e-3; % [ms] stimulus end
    
    I_parameters = [I3_mu I3_sigma I3_tau I3_onset I3_end];
    
    %% ------------------------ Noisy staircase current 
else
    % --------- Soma, input current
    I3_mustart = 0.2e-6; % [mA] mu start - initial input current
    I3_dmu = 0.05e-6; % [mA] delta mu - current step amplitude
    I3_step = 2000e-3; % [s] step duration - constant mean current during this time
    I3_sigma = 0.2e-6; % [mA] current standard deviation
    I3_zero = 0; % if true(=1), the somatic input current is equal to zero
    if I3_zero
        I3_sq = 0.*Inoise_incMu(tspan, tf, dt, I3_mustart, I3_dmu, I3_step, I3_sigma);
    else
        [I3_sq, mean_u] = Inoise_incMu(tspan, tf, dt, I3_mustart, I3_dmu, I3_step, I3_sigma);
    end
    
    % ---------- Apical, input current
    I2_mustart = 0.2e-6; % [mA] mu start - initial input current
    I2_dmu = 0.05e-6; % [mA] delta mu - current step amplitude
    I2_step = 2000e-3; % [s] step duration - constant mean current during this time
    I2_sigma = 0.09e-6; % [mA] current standard deviation
    I2_zero = 1; % if true(=1), the dendritic input current is equal to zero
    if I2_zero
        I2_sq = 0.*Inoise_incMu(tspan, tf, dt, I2_mustart, I2_dmu, I2_step, I2_sigma);
    else
        [I2_sq, mean_u] = Inoise_incMu(tspan, tf, dt, I2_mustart, I2_dmu, I2_step, I2_sigma);
    end
    
    I_parameters = [];
end

%% assigning values to Theta (parameters)
theta = [unit_number Rd_i Cai_inf I_parameters];

%% Find the steady state

if trial==1
    X0 = fsolve(@(X)FIh(0,X, theta), X0);
    save X0Ih X0
else
    load X0Ih X0
end

%% solve

if stochastic
    % Use Euler-Maruyama to integrate SDEs
    tic; [Xs, ~] = sde_euler(@(t,X)FIh(t,X,theta), @g, tspan, X0); toc
    ts = tspan;
else
    tic; [Xs, ~] = sde_euler(@(t,X)FIh(t,X,theta), @g, tspan, X0); toc
    ts = tspan;
end

%% assigning values to SV
ii=1;
Vs = Xs(:,ii:unit_number); % [mV] membrane potential at the soma PC
Vd = Xs(:,ii*unit_number+1:(ii+1)*unit_number);% ii=ii+1; % [mV] membrane potential at the tuft dendrites PC
% Ca_i = Xs(1:CSR:end,7*unit_number+1:(7+1)*unit_number) + stdInstCa.*randn(size(Xs(1:CSR:end,7),1),1);% intracellular Ca2+ concentration, it's downsampled

%% compute non state variables
Ys = zeros(numel(ts), non_state_num);
% non_state_vars = [Is, IsDx, Ib, Imbd, Itd]; point sources per neuron

for ii = 1:numel(ts)
    FIh(tspan(ii), Xs(ii,:), theta);
    Ys(ii,:) = non_state_vars;
end

% down sampling the current
Ysd = Ys(1:100:end,:);

%% Visualization
if trial<=2
    figure; plot(ts, Vs, 'k', 'LineWidth', 1); hold on; plot(ts, Vd, 'r','LineWidth', 1); xlabel('Time (s)'); ylabel('V (mV)');
    figure; plot(ts,I3_Final,'k');
    % figure; hold on; plot(ts(1:CSR:end), Ca_i); xlabel('time (s)');
    % ylabel('[Ca] (mM)');
end

%% save data

voltage.Vs = Vs;
voltage.Vd = Vd;

cd ../
file = ['GenerateFigures\FiguresData\' fileL5PCS];

if ~exist(file, 'dir') % checks if the folder already exists
    mkdir(file); 
end

save([file '\voltageTrial' num2str(trial) '.mat'], 'voltage', '-v7.3');
save([file '\currentsTrial' num2str(trial) '.mat'], 'Ysd', 'I_trans', '-v7.3');
save([file '\ICaLTrial' num2str(trial) '.mat'], 'ICaL', '-v7.3');

cd L5PCs
end

%% computing the LFPs

clearvars -except tf ntrials fileL5PCS

[lfp, lfpT] = lfp_calculation(tf, 1, ntrials, fileL5PCS);

%% Computing the CSD

iCSD_Pettersen(tf, 1, lfp, lfpT, fileL5PCS)

