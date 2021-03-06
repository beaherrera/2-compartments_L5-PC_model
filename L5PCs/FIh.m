function dXdt = FIh(t,X, theta)

%% Global variables
global non_state_vars non_state_num
global I3_sq I2_sq I3pre
global I_trans ICaL
global const_current pulse_train noisy_current
global I2_Final I3_Final dt

% created by Herrera et al., 2020

%% model by Herrera et al., 2020
jj = 1;
unit_number = theta(jj); jj=jj+1; % number of neurons;

%% assigning values to SV
ii = 1;
Vs(1:unit_number) = X(ii:unit_number); % [mV] membrane potential at the soma compartment
Vd(1:unit_number) = X(ii*unit_number+1:(ii+1)*unit_number); ii=ii+1; % [mV] membrane potential at the dendritic compartment
NamS(1:unit_number) = X(ii*unit_number+1:(ii+1)*unit_number); ii=ii+1; % Na activation gate
NahS(1:unit_number) = X(ii*unit_number+1:(ii+1)*unit_number); ii=ii+1; % Na inactivation gate
KmS(1:unit_number) = X(ii*unit_number+1:(ii+1)*unit_number); ii=ii+1; % Kdr inactivation gate
CamD(1:unit_number) = X(ii*unit_number+1:(ii+1)*unit_number); ii=ii+1; % CaL activation gate
Ca_i(1:unit_number) = X(ii*unit_number+1:(ii+1)*unit_number); ii=ii+1; % intracellular [Ca2+]
NapmD(1:unit_number) = X(ii*unit_number+1:(ii+1)*unit_number); ii=ii+1; % Nap activation gate
NaphD(1:unit_number) = X(ii*unit_number+1:(ii+1)*unit_number); ii=ii+1; % Nap inactivation gate
KsmD(1:unit_number) = X(ii*unit_number+1:(ii+1)*unit_number); ii=ii+1; % Ks activation gate
KshD(1:unit_number) = X(ii*unit_number+1:(ii+1)*unit_number); ii=ii+1; % Ks inactivation gate
Ihm(1:unit_number) = X(ii*unit_number+1:(ii+1)*unit_number);ii=ii+1; % h-current inactivation gate
Imm(1:unit_number) = X(ii*unit_number+1:(ii+1)*unit_number); % M-current activation gate

%% ------------------------------------------------------ model parameters
VshD = 8; % [mV] Shift in the dendritic channels kinetics to account for the shift in the dendritic resting potential

R_m = 50e6; %[ohm] membrane resistance soma compartment
Rd_m = 43e6; %[ohm] membrane resistance apical tuft compartment

Rd_i = theta(jj); jj=jj+1; %[ohm] transfer resistance

Cm_s = 0.26e-9; % [F] membrane capacitance soma compartment
Cm_d = 0.12e-9; % [F] membrane capacitance apical tuft compartment

%% Ionic channels soma compartment
% Na
gbar_NaS = 18e-6; % [mS] maximum conducatance of Na+ channels
E_NaS = 50; % [mV] reversal potatial for the sodium channels
[NamooS, NataumS, NahooS,NatauhS] = NaKinetics(Vs);

% K
gbar_KS = 5e-6; % [mS] maximum conducatance of K+ channels
E_KS = -85; % [mV] reversal potatial for the K+ channels
[KmooS, KtaumS] = KKinetics(Vs);

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
E_CaD = Ca_const.*log(Ca_o./Ca_i); % Nernst Equation

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
E_Ih = -45; % [mV] reversal potential for the Ih current
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
    
    % Soma
    I3 = I3_square(t, I3_sq, dt);

    % Apical dendrites
    I2 = I2_square(t, I2_sq, dt);
    
    %% Noisy Current
    
elseif noisy_current
    % ---------- Soma, input current
    I3_mu = theta(jj); jj=jj+1; % current effective amplitude
    I3_sigma = theta(jj); jj=jj+1; % [mA] noisy current standard deviation
    I3_tau = theta(jj); jj=jj+1; % [s]correlation length
    I3_onset = theta(jj); jj=jj+1;  % [s] stimulation onset
    I3_end = theta(jj); % [s] stimulation end 
    
    I3 = zeros(1,unit_number);

    % Soma
    if t==0
        I3pre = I3;
    end 
    if t>=I3_onset && t<=I3_end
        for ii=1:unit_number
            I3(ii) = Inoise(dt, I3_mu*(randn), I3_sigma*(randn), I3_tau, I3pre(ii));
            I3pre(ii) = I3(ii);
        end
    end

    I2 =0;
    
else 
    
    k = round(t./dt + 1); % index
    
    % Soma
    I3 = I3_sq(k);
    
    % Apical dendrites
    I2 = I2_sq(k);
    
end

%% Differential equations
dXdt=zeros(numel(X),1);    % creating dXdt with the number of elements of X -> No. of state variables to be solved for

%% Layer V PC
%--- currents
INa = gbar_NaS .* NamS.^3 .* NahS .* (Vs - E_NaS); % Na current
IKdr = gbar_KS.*(KmS.^4).*(Vs - E_KS); % K current
Ik_s = -INa - IKdr; % ionic currents in the soma

ICa = gbar_CaD .* (CamD.^2) .* (Vd - E_CaD); % Ca current
INap = gbar_NapD .* NapmD.^3 .* NaphD .* (Vd - E_NaD); % Nap current
IKs = gbar_KsD .* KsmD.^2.* KshD.* (Vd - E_KD); % Ks current
Ih = gbar_Ih.*Ihm.*(Vd - E_Ih); % h-current
Im = gbar_Im.*Imm.*(Vd - E_KD); % m-current
Ik_d = - ICa - INap - IKs - Ih - Im; % ionic current in the dendrites

% saving the input current applied to each compartment
indstim = round(t./dt + 1); % index
I2_Final(indstim) = I2; % dendrites
I3_Final(indstim,:) = I3;  % soma
I_trans(indstim,:) = (Vd-Vs)./Rd_i; % saving the transfer current
ICaL(indstim,:) = ICa; % saving the Ca current

%---- equations
dVsdt = ((-31.5-Vs)./R_m + (Vd-Vs)./Rd_i + Ik_s + I3)./Cm_s; % Membrane potential at the soma of layer V PC
dVddt  = ((-48.1-Vd)./Rd_m + (Vs-Vd)./Rd_i + Ik_d + I2)./Cm_d; % Membrane potential at the dendrite of layer V PC
% soma
dNamSdt = (NamooS - NamS)./NataumS; % Na activation soma
dNahSdt = (NahooS - NahS)./NatauhS; % Na inactivation soma
dKmSdt = (KmooS - KmS)./KtaumS; % K inactivation
% dendrites
dCamDdt = (CamooD - CamD)./CataumD; % Ca activation dendrite

Ad = 9302.3e-8; % [m2] area of the apical dendrites 
d = 0.1; % [um] diameter of the submembrane Ca2+ shell
gamma = 0.02; % percent of free calcium (not buffered)
factorCa = 1e4*gamma/(Ad*d);
Vd_0 = -55;
[CamD_0, ~] = CaLKineticsW(Vd_0-VshD); % Ca current at resting potential
I_Ca0 = gbar_CaD.* CamD_0.^2 .* (Vd_0 - E_CaD); % total Ca2+ current
dCa_idt = -factorCa.*((ICa-I_Ca0)./(z_Ca.*F)) - (Ca_i - Cai_inf)./tau_R; % calcium concentration in the tuft of the PC

dNapmDdt = (NapmooD-NapmD)./NaptaumD; % Nap activation dendrites
dNaphDdt = (NaphooD-NaphD)./NaptauhD; % Nap inactivation dendrites

dKsmDdt = (KsmooD-KsmD)./KstaumD; % Ks activation dendrites
dKshDdt = (KshooD-KshD)./KstauhD; % Ks inactivation dendrites

dmIhdt = (Ihmoo-Ihm)./Ihtaum; % h-current inactivation dendrites

dmImdt = (Immoo-Imm)./Imtaum; % m-current activation dendrites

%% assigning output

ii = 1;
dXdt(ii:unit_number) = dVsdt;
dXdt(ii*unit_number+1:(ii+1)*unit_number) = dVddt; ii = ii+1;
dXdt(ii*unit_number+1:(ii+1)*unit_number) = dNamSdt; ii = ii+1;
dXdt(ii*unit_number+1:(ii+1)*unit_number) = dNahSdt; ii = ii+1;
dXdt(ii*unit_number+1:(ii+1)*unit_number) = dKmSdt; ii=ii+1;
dXdt(ii*unit_number+1:(ii+1)*unit_number) = dCamDdt; ii = ii+1;
dXdt(ii*unit_number+1:(ii+1)*unit_number) = dCa_idt; ii = ii+1;
dXdt(ii*unit_number+1:(ii+1)*unit_number) = dNapmDdt; ii = ii+1;
dXdt(ii*unit_number+1:(ii+1)*unit_number) = dNaphDdt; ii = ii+1;
dXdt(ii*unit_number+1:(ii+1)*unit_number) = dKsmDdt; ii = ii+1;
dXdt(ii*unit_number+1:(ii+1)*unit_number) = dKshDdt; ii = ii+1;
dXdt(ii*unit_number+1:(ii+1)*unit_number) = dmIhdt; ii = ii+1;
dXdt(ii*unit_number+1:(ii+1)*unit_number) = dmImdt; 

%% non state variabels
% scaling factors for the returning currents, estimated from simulations
% using the detailed L5-PC model by Hay et al. 2011
% parameters names and values are in correspondance with Fig.1 from Herrera et al. 2020 JNeurosci 
alphaS_1 = 0.31682; % basal dendrites
alphaS_2 = 0.035514; % soma-AIS
alphaS_3 = 0.64767; % oblique dendrites
alphaD_1 = 0.17774; % distal trunk
alphaD_2 = 0.82226; % apical tuft
alpha_Kdr = 0.5; % scaling factor for the Kdr current

Is = (alphaS_3.*(dVsdt.*Cm_s-(-31.5-Vs)./R_m) + alpha_Kdr.*IKdr); % 3 - oblique dendrites
IsDx= (INa + alphaS_2.*(dVsdt.*Cm_s-(-31.5-Vs)./R_m)); % 2 - axon intial segment 
Ib = alphaS_1.*(dVsdt.*Cm_s-(-31.5-Vs)./R_m) + (1-alpha_Kdr).*IKdr - I3; % 1 - basal dendrites
Imbd = ICa + IKs + alphaD_1.*((dVddt.*Cm_d) - (-48.1-Vd)./Rd_m); % 4 - distal trunk
Itd = (Ih + Im + INap - I2  + alphaD_2.*((dVddt.*Cm_d) - (-48.1-Vd)./Rd_m)); % 5 - tuft dendrites 

non_state_vars = [Is, IsDx, Ib, Imbd, Itd];
non_state_num = numel(non_state_vars);

end
