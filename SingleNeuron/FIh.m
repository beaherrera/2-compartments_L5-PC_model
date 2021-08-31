function dXdt = FIh(t,X, theta)

% created by Herrera et al., 2020

%% Global variables
global non_state_vars non_state_num
global I3_sq I2_sq 
global const_current pulse_train
global I2_Final I3_Final dt

%% model by Herrera et al., 2020

%% assigning values to SV
ii = 1;
Vs = X(ii); ii=ii+1; % [mV] membrane potential at the soma compartment
Vd = X(ii); ii=ii+1; % [mV] membrane potential at the dendritic compartment
NamS = X(ii); ii=ii+1; % Na activation gate
NahS = X(ii); ii=ii+1; % Na inactivation gate
KmS = X(ii); ii=ii+1; % Kdr inactivation gate
CamD = X(ii); ii=ii+1; % CaL activation gate
Ca_i = X(ii); ii=ii+1; % intracellular [Ca2+]
NapmD = X(ii); ii=ii+1; % Nap activation gate
NaphD = X(ii); ii=ii+1; % Nap inactivation gate
KsmD = X(ii); ii=ii+1; % Ks activation gate
KshD = X(ii); ii=ii+1; % Ks inactivation gate
Ihm = X(ii); ii=ii+1; % h-current inactivation gate
Imm = X(ii); % M-current activation gate

%% ------------------------------------------------------ model parameters
VshD = 8; % [mV] Shift in the dendritic channels kinetics to account for the shift in the dendritic resting potential

R_m = 50e6; %[ohm] membrane resistance soma compartment
Rd_m = 43e6; %[ohm] membrane resistance apical tuft compartment

jj = 1;
Rd_i = theta(jj); jj=jj+1; %[ohm] transfer resistance

Cm_s = 0.26e-9; % [F] membrane capacitance soma compartment
Cm_d = 0.12e-9; % [F] membrane capacitance apical tuft compartment

%% Ionic channels soma compartment
% Na
gbar_NaS = 18e-6; % [mS] maximum conducatance of Na+ channels
E_NaS = 50; % [mV] reversal potatial for the Na+ channels
[NamooS, NataumS, NahooS,NatauhS] = NaKinetics(Vs); % gating variables and their time constant

% K
gbar_KS = 5e-6; % [mS] maximum conducatance of K+ channels 
E_KS = -85; % [mV] reversal potatial for the K+ channels
[KmooS, KtaumS] = KKinetics(Vs); % gating variables and their time constant

%% ionic channels dendritic compartment
% Calcium L-type - CaL 
gbar_CaD = 3.85e-6; % [mS] maximum conducatance for Ca2+ L type channels
% Equilibrium Potential Ca
R = 8.314e3; % [mJK-1mol-1] gas constant
z_Ca = 2; % valence of Ca2+
F = 9.648e4; % [Cmol-1] Faraday's constant
temp = 310.15; % [K] temperature
tau_R = 80e-3; % [s] Ca decay time constant
Ca_o = 2; % [mM] extracellular calcium concentration
Cai_inf = theta(jj); jj=jj+1; % [mM] intracellular equilibrium calcium concentration
Ca_const = (R*temp)/(z_Ca*F);
E_CaD = Ca_const.*log(Ca_o./Ca_i); %Equilibrium Potential Ca

[CamooD, CataumD] = CaLKineticsW(Vd-VshD); % gating variables and their time constant

% persistent sodium - Nap
gbar_NapD = 0.022e-6; % [S] maximum conducatance for persistent Na+ channels
E_NaD = 50; % [mV] reversal potatial for the sodium channels
[NapmooD, NaptaumD, NaphooD, NaptauhD] = NapKinetics(Vd-VshD); % gating variables and their time constant

% slow inactivating potasium - Ks
gbar_KsD = 28e-6; % [S] maximum conducatance for slow K+ channels
E_KD = -85; % [mV] reversal potatial for the sodium channels
[KsmooD, KstaumD, KshooD, KstauhD] = KslowKinetics(Vd-VshD); % gating variables and their time constant

% Ih current
E_Ih = -45; % [mV] reversal potentail used for the Ih current
gbar_Ih = 0.865e-6; % [S] maximum conductance of non-specific cation channels
[Ihmoo,Ihtaum]=IhKinetics(Vd-VshD); % gating variables and their time constant

% Im current
gbar_Im = 1e-6; % [S] maximum conducatance
[Immoo,Imtaum]=ImKinetics(Vd-VshD); % gating variables and their time constant

%% --------------------------------- stimulation protocol
%% constant current

if const_current
    % CI
    I2 = 0; % [mA] amplitude of the dendritic input I2
    I3 = 0; % [mA] amplitude of the somatic input I3
    
    %% Soma
    I3_amp = theta(jj); jj=jj+1; % [mA] amplitude of the somatic input I3
    I3_stim_onset = theta(jj); jj=jj+1; % [s] stimulus onset
    I3_stim_end = theta(jj); jj=jj+1; % [s] stimulus end
    
    if t>=I3_stim_onset && t<=I3_stim_end
        I3 = I3_amp; % [mA]
    end
    
    %% Apical dendrites
    I2_amp = theta(jj); jj=jj+1; % [mA] amplitude of the dendritic input I2
    I2_stim_onset = theta(jj); jj=jj+1; % [s] stimulus onset
    I2_stim_end = theta(jj); % [s] stimulus end
    
    if I2_stim_end==0 % if I2_stim_end is equal to zero, a EPSP-like stimuli is generated
        if t>=I2_stim_onset
            I2 = I2_amp.*doubleExpFunct(t - I2_stim_onset, 2e-3, 8e-3);
        end
    else % otherwise, a current step pulse is generated
        if t>=I2_stim_onset && t<=I2_stim_end
            I2 = I2_amp; % [mA]
        end
        
    end
    
    %% Train of pulses
    
elseif pulse_train
    
    k = round(t./dt + 1); % index
    
    % Soma
    I3 = I3_sq(k);
    
    % Apical dendrites
    I2 = I2_sq(k);
    
     %% Noisy staircase current 
else 
    
    k = round(t./dt + 1); % index
    
    % Soma
    I3 = I3_sq(k);
    
    % Apical dendrites
    I2 = I2_sq(k);
    
end

%% Differential equations

dXdt=zeros(numel(X),1); % creating dXdt with the number of elements of X -> No. of state variables to be solved for

%% Layer V PC
%--- currents
INa = gbar_NaS .* NamS.^3 .* NahS .* (Vs - E_NaS); % Na current
IKdr = gbar_KS.*(KmS.^4).*(Vs - E_KS); % K current
Ik_s = -INa - IKdr; % ionic currents in the soma

ICa = gbar_CaD .* (CamD.^2) .* (Vd - E_CaD);  % CaL current
INap = gbar_NapD .* NapmD.^3 .* NaphD .* (Vd - E_NaD); % Nap current
IKs = gbar_KsD .* KsmD.^2.* KshD.* (Vd - E_KD); % Ks current
Ih = gbar_Ih.*Ihm.*(Vd - E_Ih); % h-current
Im = gbar_Im.*Imm.*(Vd - E_KD); % m-current
Ik_d = - ICa - INap - IKs - Ih - Im; % ionic current in the dendrites

% saving the input current applied to each compartment
indstim = round(t./dt + 1);
I2_Final(indstim) = I2;
I3_Final(indstim) = I3;

%---- equations
dVsdt = ((-31.5-Vs)./R_m + (Vd-Vs)./Rd_i + Ik_s + I3)./Cm_s; % Membrane potential at the soma
dVddt  = ((-48.1-Vd)./Rd_m + (Vs-Vd)./Rd_i + Ik_d + I2)./Cm_d; % Membrane potential at the dendrite

 % soma
dNamSdt = (NamooS - NamS)./NataumS; % INa activation soma
dNahSdt = (NahooS - NahS)./NatauhS; % INa inactivation soma
dKmSdt = (KmooS - KmS)./KtaumS; % K inactivation

 % dendrites
dCamDdt = (CamooD - CamD)./CataumD; % ICa activation dendrite

Ad = 9302.3e-8; % [m2] area of the apical dendrites 
d = 0.1; % [um] diameter of the submembrane Ca2+ shell
gamma = 0.02; % percent of free calcium (not buffered)
factorCa = 1e4*gamma/(Ad*d);
Vd_0 = -55;
[CamD_0, ~] = CaLKineticsW(Vd_0-VshD);
I_Ca0 = gbar_CaD.* CamD_0.^2 .* (Vd_0 - E_CaD); % Ca current at resting potential
dCa_idt = -factorCa.*((ICa-I_Ca0)./(z_Ca.*F)) - (Ca_i - Cai_inf)./tau_R; % calcium concentration in the dendritic compartment

dNapmDdt = (NapmooD-NapmD)./NaptaumD; % Nap activation dendrites
dNaphDdt = (NaphooD-NaphD)./NaptauhD; % Nap inactivation dendrites

dKsmDdt = (KsmooD-KsmD)./KstaumD; % Ks activation dendrites
dKshDdt = (KshooD-KshD)./KstauhD; % Ks inactivation dendrites

dmIhdt = (Ihmoo-Ihm)./Ihtaum; % h-current inactivation dendrites

dmImdt = (Immoo-Imm)./Imtaum; % m-current activation dendrites

%% assigning output
ii = 1;
dXdt(ii) = dVsdt; ii = ii+1;
dXdt(ii) = dVddt; ii = ii+1;
dXdt(ii) = dNamSdt; ii = ii+1;
dXdt(ii) = dNahSdt; ii = ii+1;
dXdt(ii) = dKmSdt; ii=ii+1;
dXdt(ii) = dCamDdt; ii = ii+1;
dXdt(ii) = dCa_idt; ii = ii+1;
dXdt(ii) = dNapmDdt; ii = ii+1;
dXdt(ii) = dNaphDdt; ii = ii+1;
dXdt(ii) = dKsmDdt; ii = ii+1;
dXdt(ii) = dKshDdt; ii = ii+1;
dXdt(ii) = dmIhdt; ii = ii+1;
dXdt(ii) = dmImdt; 

%% non state variabels

non_state_vars = [INa, IKdr, ICa, INap, IKs, Ih, Im];
non_state_num = numel(non_state_vars);

end
