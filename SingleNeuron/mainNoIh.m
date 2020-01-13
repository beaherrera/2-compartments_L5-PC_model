%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Layer 5 Pyramidal cell mininal model
%
% Compartments: basal-dendritic/somatic compartment and
% apical-dendrite/trunk compartment
% - Soma/basal-dendrites compartment includes typical Na+ and K+ conductances. 
% - Trunk/apical-dendrites compartments includes: persistent Na+ (Nap), 
%   hyperpolarization-activated cation (Ih), slow inactivation K+ (Ks), 
%   muscarinic K+ (IM) and Ca2+ L-type currents (CaL). 
%
% In this program, the h-current was blocked.
%
% model by Herrera et al., 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%% Global variables
global SigX const_current pulse_train
global I3_sq I2_sq I2_Final I3_Final
global dt
% global non_state_vars non_state_num % uncomment if you want to comput the
% ionic currents. Uncomment the part of the non-state variables calculation
% as well.

%% simulation time
t0 = 0; %[s] start time
tf = 110*1e-3;  %[s] end time 
dt = 0.000001;  %[s] time increment
CSR=round(0.001./dt); % downsampling the tspan to export [Ca2+]i  

tspan = t0:dt:tf;   %[s] time span
tspan = tspan(1:(numel(tspan)-1));

%% stochastic or deterministic
stochastic = true;  % if true system is stochastic - if false system is deterministic

%%  ----------------------------------------------- parameters

Rd_i = 65e6; % [Ohm] transfer resistance between the Soma/basal-dendrites and Trunk/apical-dendrites compartments
plotresults = 1; % if plotresults=1; figures with the voltages at both compartments 
%                   and the [Ca2+]i traces will be generated after solving the model
f = 105; % frequency of the train of pulses when pulse_train =1

%% --------------------------------------------- Stimulus
% pre-allocating memory
% these variables store the dendritic and somatic stimulus for later
% visualization
I2_Final = zeros(size(tspan)); % somatic input current 
I3_Final = I2_Final; % dendritic input current

% stimulus type
const_current = 1; % if true(=1): soma stimuli is a current step pulse, and the
%                    dendritic stimuli is either a current step pulse or a double exponential
pulse_train = 0; % if true(=1): the somatic/dendritic input is a train of square current pulses
% When both of these variables are false(=0) the somatic/dendritic input
% current is a Noisy staircase current. This type of current was used to
% generate the f-I curve for either somatic or dendritic stimulations.

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
    % spike generation. We applied current pulses with the following
    % frequencies: [30:10:90, 105, 110:10:170] Hz.
    % The critical frequency when Ih is blocked was 105 Hz.
    
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

%% ----------------------------------------------- initial conditions (ICs)
% ICs are set here. However after defining the parameters in the next
% section, the steady state of the system is found using fsolve and the
% results will be used as the actual ICs.

Vs_0 = -65; % somatic voltage
Vd_0 = -65; % dendritic voltage
[NamS_0, ~, NahS_0, ~] = NaKinetics(Vs_0); % Na gating variables
[KmS_0,~]=KKinetics(Vs_0); % K gating variables
[CamD_0, ~] = CaLKineticsW(Vd_0); % CaL gating variables
Cai_inf = 80e-6; % [mM] intracellular equilibrium calcium concentration at rest
Cai_0 = Cai_inf; % [mM] 
[NapmD_0, ~, NaphD_0, ~] = NapKinetics(Vd_0); % Nap gating variables
[KsmD_0, ~, KshD_0, ~] = KslowKinetics(Vd_0); % Ks gating variables
[Imm_0, ~]=ImKinetics(Vd_0); % Im gating variables

X0 = [Vs_0, Vd_0, NamS_0, NahS_0, KmS_0, CamD_0, Cai_0, NapmD_0, NaphD_0, KsmD_0, KshD_0, Imm_0]';

state_var_num = numel(X0);  % number of state variabels

%% assigning values to Theta (parameters)
theta = [Rd_i Cai_inf I_parameters];

%% Find the steady state and set it as initial conditions

X0 = fsolve(@(X)FnoIh(0,X, theta), X0);

%% Variance of Wiener Process
% We considered only the intracellular calcium concentration as an
% stochastic process.
g_Vs = 0;
g_Vd = 0;
g_NamS = 0;
g_NahS = 0;
g_KmS = 0;
g_CamD = 0;
g_Cai = 0.000000001;
g_NapmD = 0;
g_NaphD = 0;
g_KsmD = 0;
g_KshD = 0;
g_ImmD = 0;

stdInstCa = 0.00001; % intrumental noise of the Ca2+ traces

if stochastic
    SigX = [g_Vs, g_Vd, g_NamS, g_NahS, g_KmS, g_CamD, g_Cai, g_NapmD, g_NaphD, g_KsmD, g_KshD, g_ImmD]'./sqrt(dt);
else
    SigX = zeros(state_var_num,1);
end

%% solve

if stochastic
   
    % Use Euler-Maruyama to integrate SDEs
    tic; [Xs, ~] = sde_euler(@(t,X)FnoIh(t,X,theta), @g, tspan, X0); toc
    ts = tspan;
else
    tic; [ts,Xs] = ode452(@(t,X)FnoIh(t,X,theta), [t0, tf], X0); toc
end

%% assigning values to state variables
ii=1;
Vs = Xs(:,ii); ii= ii+1; % [mV] membrane potential at the soma PC
Vd = Xs(:,ii); % [mV] membrane potential at the tuft dendrites PC
Ca_i = Xs(1:CSR:end,7) + stdInstCa.*randn(size(Xs(1:CSR:end,7),1),1); % intracellular Ca2+ concentration, it's downsampled

%% Plotting the results

if plotresults
    figure;
    plot(ts, Vs, 'k', 'LineWidth', 1);
    hold on; plot(ts, Vd, 'r', 'LineWidth', 1);
    xlabel('Time (s)'); ylabel('V (mV)');
    plot(ts,5.*I3_Final./max(abs(I3_Final))-95,'k');
    plot(ts,5.*I2_Final./max(abs(I2_Final))-90,'r');
    
    figure; hold on; plot(ts(1:CSR:end), Ca_i); xlabel('time (ms)'); ylabel('[Ca] (mM)');
end

%% Compute I-F curve

if ~const_current && ~pulse_train
    % This computes the f-I curve for dendritic or somatic stimulations
    % when a ramp is used as stimuli
    
    if nnz(I3_sq) == length(I3_sq) 
        current_freqRelation(Vs, I3_step, tf, tspan, I3_mustart, I3_dmu)
    else
        current_freqRelation(Vs, I2_step, tf, tspan, I2_mustart, I2_dmu)
    end
    
end


%% compute non state variables
% % the ionic currents are the non-state variables
% Ys = zeros(numel(ts), non_state_num);
% % non_state_vars = [INa, IKdr, ICa, INap, IKs, Im];
% 
% for ii = 1:numel(ts)
%     FnoIh(tspan(ii), Xs(ii,:), theta);
%     Ys(ii,:) = non_state_vars;
% end

%% save data

% Uncomment the code below if you want to do a loop for all the frequencies and
% store the results in a single file. The code concatenates the voltages, the [Ca2+]
% and the stimulation current. Start the loop before declering the stimulus type 
% and ended it after the condition below ends
% if pulse_train
%     if f==30
%         Vsf = Vs;
%         Vdf = Vd;
%         I3f = I3_Final;
%         Ca_if = Ca_i;
%     else
%         Vsf = cat(2, Vsf, Vs);
%         Vdf = cat(2, Vdf, Vd);
%         I3f = cat(2, I3f, I3_Final);
%         Ca_if = cat(2, Ca_if, Ca_i);
%     end
% end

cd ../
cd GenerateFigures
file = 'FiguresData';
if ~exist(file, 'dir') % checks if the folder already exists
    mkdir(file);  % creates a folder named 'file'
end

if const_current % save the simulated data for the BAC-firing paradigm
    fileB = [file '\BACfiring'];
    if ~exist(fileB, 'dir') % checks if the folder already exists
        mkdir(fileB); % creates a folder named 'fileB'
    end
    save([fileB '\voltagenoIh' num2str(paradigm) '.mat'], 'Vs', 'Vd');
    save([fileB '\stimulusnoIh' num2str(paradigm) '.mat'], 'I3_Final', 'I2_Final');
    
elseif pulse_train % save the simulated data for the critical frequency study
    fileCF = [file '\CF'];
    if ~exist(fileCF, 'dir') % checks if the folder already exists
        mkdir(fileCF); % creates a folder named 'fileCF'
    end
    save([fileCF '\voltagenoIhFreq.mat'], 'Vsf', 'Vdf', '-v7.3');
    save([fileCF '\stimulusnoIhFreq.mat'], 'I3f', 'I2_Final', '-v7.3');
    save([fileCF '\calciumnoIhFreq.mat'], 'Ca_if');
    
else % save the simulated data for a ramp input current
    fileR = [file '\CurrentRamp'];
    if ~exist(fileR, 'dir') % checks if the folder already exists
        mkdir(fileR); % creates a folder named 'fileR'
    end
    if nnz(I2_sq) ~= length(I2_sq) % True for somatic stimulation
        save([fileR '\voltagenoIhRamp.mat'], 'VsR', 'VdR', '-v7.3');
        save([fileR '\stimulusnoIhRamp.mat'], 'I3_Final', 'mean_u', '-v7.3');
    else % True for dendritic stimulation
        save([fileR '\voltagenoIhRampD.mat'], 'VsR', 'VdR', '-v7.3');
        save([fileR '\stimulusnoIhRampD.mat'], 'I2_Final', 'mean_u', '-v7.3');
    end
end

