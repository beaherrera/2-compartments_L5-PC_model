%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% fig_CritFreq
% Authors: B. Herrera et al. 2020
% This program generates the Figure 4. Variation of dendritic Ca2+ spikes 
% occurrence relative to somatic stimulation frequency. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
clc

%% simulation time
t0 = 0; %[s] start time
tf = 50*1e-3;  %[s] end time 
dt = 0.000001;  %[s] time increment

tspan = t0:dt:tf;   %[s] time span
tspan = tspan(1:(numel(tspan)-1));

% Ca2+ dynamics time span
CSR=round(0.001./dt);
tspanCa1 = t0:dt:110e-3;
tspanCa = tspanCa1(1:CSR:end-1);

%% -- loading simulated data
% blocking Ih
load('FiguresData\CF\voltagenoIhFreq.mat', 'Vdf');
VdfnoIh = Vdf;
clear Vdf
% Ih present
load('FiguresData\CF\voltageIhFreq.mat', 'Vsf', 'Vdf');
load('FiguresData\CF\calciumIhFreq.mat', 'Ca_if');

%% computing the integral of the dendritic voltage traces

VdIntIh = trapz(tspan(1:length(tspan)-1000), Vdf(1:length(tspan)-1000,:)+55); % Ih present
VdIntnoIh = trapz(tspan(1:length(tspan)-1000), VdfnoIh(1:length(tspan)-1000,:)+65); % Ih blocked

% Note: the respective dendritic resting potential were added in each case
% to compute the integral relative to 0mV to obtain positive values.

%% visualization

figure('units','inch','position',[0,0,6,6]);  % defining the figure's dimensions
font = 10; % font size

fIh=[30:10:140 149 160 170]; % stimulation frequencies when Ih is present
fnoIh=[30:10:90, 105, 110:10:170]; % stimulation frequencies when Ih is blocked

CFIh = 149; % critical frequency when Ih is present
CFnoIh = 105; % critical frequency when Ih is blocked

% =========================== A
%====== below critical frequency
% --- Stimulus
pos1 = [0.05 0.9 0.43 0.05]; % position of the plot
subplot('Position',pos1)
% input current
len = length(tspan);
I3_amp = 15.*1e-6; % [mA] amplitude of the somatic input current, I3
widthAct_I3 = 2e-3; % [s] width of the pulses
Freq_I3 = 100; % [Hz] pulses frequency
I3_stim_onset = 1e-3; % [s] stimulus onset
I3_stim_end = 110e-3; % [s] stimulus end
PWindowI3 = round([I3_stim_onset,I3_stim_end]./dt);
width_I3 = round(widthAct_I3./dt); % width of the square pulse in time pts
IP_I3 = round(1./(Freq_I3.*dt) - width_I3); % inter-pulse interval
NC_I3 = round((tf./dt)./(width_I3 + IP_I3)); % number of pulses
I3_sq = square_pulse(I3_amp, width_I3, IP_I3, NC_I3, len, PWindowI3); % generates the square pulses
% plotting the input current
plot(tspan(1:length(tspan)-2000), I3_sq(1:(length(tspan)-2000)).*1e6, '-k', 'LineWidth', 1); 
hold on
plot([0; 0]+45e-3, [2; 12], '-k', 'LineWidth', 1.5) % voltage scale
hold off
text(0+0.1e-3+45e-3, mean([2 12]), ' 10 nA', 'HorizontalAlignment','left', 'fontsize', font) % voltage scale label
text(6e-3, 7, 'I_{stim}', 'HorizontalAlignment','right', 'fontsize', font,'Color','k')
text(0, 20, 'A', 'HorizontalAlignment','right', 'fontsize', font+4,'fontweight','bold') % letter of the figure
text(mean([0 tspan(end)]), 25, '100 Hz', 'HorizontalAlignment','center', 'fontsize', font)
set(gca, 'Visible', 'off')

% --- Somatic voltage
pos2 = [0.05 0.73 0.43 0.15]; % position of the plot
subplot('Position',pos2)
% plotting the voltage
plot(tspan(1:length(tspan)-1800), Vsf(1:(length(tspan)-1800),fIh==100), '-k', 'LineWidth', 1); 
hold on
plot([0; 0]+45e-3, [-30; 10], '-k', 'LineWidth', 1.5) % voltage scale
hold off
text(0+1e-3+45e-3, mean([-30; 10]), '50 mV', 'HorizontalAlignment','left', 'fontsize', font) % voltage scale label
text(3e-3, -40, 'V_s', 'HorizontalAlignment','right', 'fontsize', font,'Color','k')
set(gca, 'Visible', 'off')

% --- dendritic voltage
pos3 = [0.05 0.58 0.43 0.15]; % position of the plot
subplot('Position',pos3)
plot(tspan(1:length(tspan)-1000), Vdf(1:(length(tspan)-1000),fIh==100), '-r', 'LineWidth', 1); % dendritic voltage at 100 Hz
hold on
plot(tspan(length(tspan)-1000)-[0; 10e-3], [-65; -65], '-k', 'LineWidth', 1.5);  % time scale
hold on
% selecting the zones that correspond to the intracellular Ca2+ traces
% plotted in panel B
plot([0; 0]+45e-3, [-20; 0], '-k', 'LineWidth', 1.5) % current scale
hold on
plot([0.0084 0.0084], [min(Vdf(1:(length(tspan)-1000),fIh==CFIh)) max(Vdf(1:(length(tspan)-1000),fIh==CFIh))], '-b', 'LineWidth', 1);
hold on
plot([0.0185 0.0185], [min(Vdf(1:(length(tspan)-1000),fIh==CFIh)) max(Vdf(1:(length(tspan)-1000),fIh==CFIh))], '-b', 'LineWidth', 1);
hold on
plot([0.0279 0.0279], [min(Vdf(1:(length(tspan)-1000),fIh==CFIh)) max(Vdf(1:(length(tspan)-1000),fIh==CFIh))], '-b', 'LineWidth', 1);
hold on
plot([0.03796 0.03796], [min(Vdf(1:(length(tspan)-1000),fIh==CFIh)) max(Vdf(1:(length(tspan)-1000),fIh==CFIh))], '-b', 'LineWidth', 1);
hold off
text(tspan(length(tspan)-1000)-mean([0 10e-3]), -80, '10 ms', 'HorizontalAlignment','center', 'fontsize', font) % time scale label
text(0+1e-3+45e-3, mean([-20; 0]), '20 mV', 'HorizontalAlignment','left', 'fontsize', font) % current scale label
text(3e-3, -30, 'V_d', 'HorizontalAlignment','right', 'fontsize', font,'Color','red')
ylim([min(Vdf(1:(length(tspan)-1000),fIh==CFIh)) max(Vdf(1:(length(tspan)-1000),fIh==CFIh))])
set(gca, 'Visible', 'off')

% linking the dendritic voltage with the corresponding intracellular Ca2+ concentration in panel B-left 
h1 = annotation('line',[0.08 0.123],[0.55 0.58]);
h1.LineStyle = '--';
h1.Color = 'b';
h1.LineWidth = 1;
h2 = annotation('line',[0.12 0.21],[0.55 0.58]);
h2.LineStyle = '--';
h2.LineWidth = 1;
h2.Color = 'b';
h3 = annotation('line',[0.15 0.29],[0.55 0.58]);
h3.LineStyle = '--';
h3.LineWidth = 1;
h3.Color = 'b';
h4 = annotation('line',[0.19 0.37],[0.55 0.58]);
h4.LineStyle = '--';
h4.LineWidth = 1;
h4.Color = 'b';

% =========================== B
% --- Intracellular Calcium Concentration
pos4 = [0.05 0.4 0.43 0.15]; % position of the plot
subplot('Position',pos4)
plot(tspanCa, Ca_if(1:length(tspanCa),fIh==100), '-k', 'LineWidth', 1); % [Ca2+]i
ylim([3e-5 max(Ca_if(1:length(tspanCa),fIh==CFIh))])
hold on
plot(tspanCa(end)-[0; 10e-3], [3.5e-5; 3.5e-5], '-k', 'LineWidth', 1.5);  % time scale
hold on
plot([tspanCa(end); tspanCa(end)], [1e-4; 2e-4], '-k', 'LineWidth', 1.5) % current scale
hold on
% selecting the zones that correspond to the voltage intervals 
% plotted in panel A bottom-left
plot([0.0084 0.0084], [min(Ca_if(1:length(tspanCa),fIh==100)) max(Ca_if(1:length(tspanCa),fIh==CFIh))], '-b', 'LineWidth', 1);
hold on
plot([0.0185 0.0185], [min(Ca_if(1:length(tspanCa),fIh==100)) max(Ca_if(1:length(tspanCa),fIh==CFIh))], '-b', 'LineWidth', 1);
hold on
plot([0.0279 0.0279], [min(Ca_if(1:length(tspanCa),fIh==100)) max(Ca_if(1:length(tspanCa),fIh==CFIh))], '-b', 'LineWidth', 1);
hold on
plot([0.03796 0.03796], [min(Ca_if(1:length(tspanCa),fIh==100)) max(Ca_if(1:length(tspanCa),fIh==CFIh))], '-b', 'LineWidth', 1);
hold off
text(mean(tspanCa(end)-[0 10e-3]), 0.5e-5, '10 ms', 'HorizontalAlignment','center', 'fontsize', font) % time scale label
text(tspanCa(end)+2e-3, mean([1e-4; 2e-4]), '100 nM', 'HorizontalAlignment','left', 'fontsize', font) % current scale label
text(0, max(Ca_if(:,fIh==CFIh)), 'B', 'HorizontalAlignment','right', 'fontsize', font+4,'fontweight','bold') % letter of the figure
text(5e-3, 15e-5, '[Ca]_{i}', 'HorizontalAlignment','right', 'fontsize', font,'Color','k')
set(gca, 'Visible', 'off')

% =========================== A
%====== at critical frequency
% --- Stimulus
pos1 = [0.53 0.9 0.43 0.05]; % position of the plot
subplot('Position',pos1)
% input current
len = length(tspan);
I3_amp = 15.*1e-6;% [mA] amplitude of the somatic input current, I3
widthAct_I3 = 2e-3; % [s] width of the pulses
Freq_I3 = CFIh; % [Hz] pulses frequency
I3_stim_onset = 1e-3; % [s] stimulus onset
I3_stim_end = 110e-3; % [s] stimulus end
PWindowI3=round([I3_stim_onset,I3_stim_end]./dt);
width_I3 = round(widthAct_I3./dt); % width of the square pulse in time pts
IP_I3 = round(1./(Freq_I3.*dt) - width_I3); % inter-pulse interval
NC_I3 = round((tf./dt)./(width_I3 + IP_I3)); % number of pulses
I3_sq = square_pulse(I3_amp, width_I3, IP_I3, NC_I3, len, PWindowI3); % generates the square pulses
% plotting the input current
plot(tspan(1:length(tspan)-2000), I3_sq(1:(length(tspan)-2000)).*1e6, '-k', 'LineWidth', 1); 
text(mean([0 tspan(end)]), 25, [num2str(CFIh) ' Hz'], 'HorizontalAlignment','center', 'fontsize', font)
set(gca, 'Visible', 'off')

% --- Somatic voltage
pos2 = [0.53 0.73 0.43 0.15]; % position of the plot
subplot('Position',pos2)
% plotting the voltage
plot(tspan(1:length(tspan)-1800), Vsf(1:(length(tspan)-1800),fIh==CFIh), '-k', 'LineWidth', 1); 
set(gca, 'Visible', 'off')

% --- dendritic voltage
pos3 = [0.53 0.57 0.43 0.15]; % position of the plot
subplot('Position',pos3)
plot(tspan(1:length(tspan)-1000), Vdf(1:(length(tspan)-1000),fIh==CFIh), '-r', 'LineWidth', 1); % voltage
hold on
% selecting the zones that correspond to the intracellular Ca2+ traces
% plotted in panel B-right
plot([0.0055 0.0055], [min(Vdf(1:(length(tspan)-1000),fIh==CFIh)) max(Vdf(1:(length(tspan)-1000),fIh==CFIh))], '-b', 'LineWidth', 1);
hold on
plot([0.0128 0.0128], [min(Vdf(1:(length(tspan)-1000),fIh==CFIh)) max(Vdf(1:(length(tspan)-1000),fIh==CFIh))], '-b', 'LineWidth', 1);
hold on
plot([0.019 0.019], [min(Vdf(1:(length(tspan)-1000),fIh==CFIh)) max(Vdf(1:(length(tspan)-1000),fIh==CFIh))], '-b', 'LineWidth', 1);
hold on
plot([0.045 0.045], [min(Vdf(1:(length(tspan)-1000),fIh==CFIh)) max(Vdf(1:(length(tspan)-1000),fIh==CFIh))], '-b', 'LineWidth', 1);
hold off
set(gca, 'Visible', 'off')

% linking the dendritic voltage with the corresponding intracellular Ca2+ concentration in panel B 
h1 = annotation('line',[0.55 0.57],[0.55 0.58]);
h1.LineStyle = '--';
h1.Color = 'b';
h1.LineWidth = 1;
h2 = annotation('line',[0.578 0.64],[0.55 0.58]);
h2.LineStyle = '--';
h2.LineWidth = 1;
h2.Color = 'b';
h3 = annotation('line',[0.595 0.685],[0.55 0.58]);
h3.LineStyle = '--';
h3.LineWidth = 1;
h3.Color = 'b';
h4 = annotation('line',[0.69 0.915],[0.55 0.58]);
h4.LineStyle = '--';
h4.LineWidth = 1;
h4.Color = 'b';

% =========================== B
% --- Intracellular Calcium Concentration
pos4 = [0.53 0.4 0.43 0.15]; % position of the plot
subplot('Position',pos4)
plot(tspanCa, Ca_if(1:length(tspanCa),fIh==CFIh), '-k', 'LineWidth', 1); 
hold on
% selecting the zones that correspond to the voltage intervals 
% plotted in panel A bottom-right
plot([0.0055 0.0055], [min(Ca_if(1:length(tspanCa),fIh==CFIh)) max(Ca_if(1:length(tspanCa),fIh==CFIh))], '-b', 'LineWidth', 1);
hold on
plot([0.0128 0.0128], [min(Ca_if(1:length(tspanCa),fIh==CFIh)) max(Ca_if(1:length(tspanCa),fIh==CFIh))], '-b', 'LineWidth', 1);
hold on
plot([0.019 0.019], [min(Ca_if(1:length(tspanCa),fIh==CFIh)) max(Ca_if(1:length(tspanCa),fIh==CFIh))], '-b', 'LineWidth', 1);
hold on
plot([0.045 0.045], [min(Ca_if(1:length(tspanCa),fIh==CFIh)) max(Ca_if(1:length(tspanCa),fIh==CFIh))], '-b', 'LineWidth', 1);
hold off
set(gca, 'Visible', 'off')
ylim([3e-5 max(Ca_if(1:length(tspanCa),fIh==CFIh))])

% =========================== C
% --- Integral
pos4 = [0.1 0.08 0.35 0.25]; % position of the plot
subplot('Position',pos4)
pIh = plot(fIh, VdIntIh, '.k', 'markersize', 18); % Ih present
hold on
plot(fIh, VdIntIh, '-k', 'markersize', 0.5); % Ih present
hold on;
pnoIh = plot(fnoIh, VdIntnoIh, 'sb', 'MarkerFaceColor', 'b', 'markersize', 4); % Ih blocked
hold on
plot(fnoIh, VdIntnoIh, '-b', 'markersize', 0.5); % Ih blocked
legend([pIh pnoIh], {'Ih Present', 'Ih Blocked'},...
    'Location','northeast', 'fontsize', font-2,'linewidth',1)
legend('boxoff')
ylabel('Integral (mV*s)');
xlabel('Frequency (Hz)');
text(-30, 0.45, 'C', 'HorizontalAlignment','right', 'fontsize', font+4,'fontweight','bold') % letter of the figure
x1 = [0.235 0.265];
y2 = [0.29 0.32];
annotation('textarrow',x1,y2,'String',[num2str(CFnoIh) ' Hz'])
x1 = [0.332 0.356];
y2 = [0.15 0.19];
annotation('textarrow',x1,y2,'String',[num2str(CFIh) ' Hz'])
box 'off'
set(gca,'linewidth',1.5,'fontsize',font,'fontweight','bold')

% =========================== D
% --- Experimental data (Berger et al. 2003)
CF_Ih = [165, 138, 137, 125, 127, 110, 109, 105, 100, 80, 80]; % Ih present
CF_noIh = [120, 90, 95, 90, 80, 85, 60, 70, 65, 50, 60]; % Ih blocked

% mean CF
mean_CF_Ih = mean(CF_Ih); % Ih present
mean_CF_noIh = mean(CF_noIh); % Ih present

pos4 = [0.6 0.078 0.35 0.25]; % position of the plot
subplot('Position',pos4)
% bar plot
Y = [mean_CF_Ih mean_CF_noIh];
bar(Y,0.5,'FaceColor',[0.85 0.85 0.85])
hold on
% line plots
for ii=1:length(CF_Ih)
    if ii==1
        pExp = plot([1 2], [CF_Ih(ii) CF_noIh(ii)], 'o-k', 'markersize', 6 ,'linewidth',1);
    else
        plot([1 2], [CF_Ih(ii) CF_noIh(ii)], 'o-k', 'markersize', 6 ,'linewidth',1)
    end
end
hold on
pM = plot([1 2], [138 105], 'd-b','MarkerFaceColor', 'b', 'markersize', 7,'linewidth',1); % model predictions
hold off
xticks([1 2])
xlim([0.5 2.5])
tickslab = {'Ih Present','Ih Blocked'};
xticklabels(tickslab)
box 'off'
legend([pExp pM],{'Experiment', 'Model'},...
    'Location','northeast', 'fontsize', font-2,'linewidth',1)
set(gca,'linewidth',1.5,'fontsize',font,'fontweight','bold')
ylabel('Critical Frequency (Hz)');
text(-0.1, 220, 'D', 'HorizontalAlignment','right', 'fontsize', font+4,'fontweight','bold') % letter of the figure

%% -- saving the figure
file = 'Figures';
if ~exist(file, 'dir') % checks if the folder already exists
    mkdir(file);  % creates a folder named 'file'
end

print(gcf, '-dtiff', '-r1000', 'Figures\fig_CritFreq.tiff');

