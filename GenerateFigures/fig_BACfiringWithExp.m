%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% fig_BACfiringWithExp
% Authors: B. Herrera et al. 2020
% This program generates the Figure 3. Back-propagating AP activated Ca2+ 
% spike firing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
clc

%% simulation time
t0 = 0; %[s] start time
tf = 100*1e-3;  %[s] end time 
dt = 0.000001;  %[s] time increment

tspan = t0:dt:tf;   %[s] time span
tspan = 1e3.*tspan(1:(numel(tspan)-1)); % converting from s -> ms

%% -- loading simulated data
 %   paradigm=1; subthreshold dendritic stimulation only
load('FiguresData\BACfiring\voltageIh1.mat', 'Vs', 'Vd');
load('FiguresData\BACfiring\stimulusIh1.mat', 'I3_Final', 'I2_Final');
Vs1 = Vs(1:length(tspan))+65; % voltage at the soma
Vd1 = Vd(1:length(tspan))+55; % voltage at the dendrites
Is1 = I3_Final(1:length(tspan)); % somatic input 
Id1 = I2_Final(1:length(tspan)); % dendritic input
%   paradigm=2; suprathreshold somatic stimulation only
load('FiguresData\BACfiring\voltageIh2.mat', 'Vs', 'Vd');
load('FiguresData\BACfiring\stimulusIh2.mat', 'I3_Final', 'I2_Final');
Vs2 = Vs(1:length(tspan))+65; % voltage at the soma
Vd2 = Vd(1:length(tspan))+55; % voltage at the dendrites
Is2 = I3_Final(1:length(tspan)); % somatic input 
Id2 = I2_Final(1:length(tspan)); % dendritic input
%   paradigm=3; combine subthreshold dendritic and suprathreshold
%               somatic stimulations
load('FiguresData\BACfiring\voltageIh3.mat', 'Vs', 'Vd');
load('FiguresData\BACfiring\stimulusIh3.mat', 'I3_Final', 'I2_Final');
Vs3 = Vs(1:length(tspan))+65; % voltage at the soma
Vd3 = Vd(1:length(tspan))+55; % voltage at the dendrites
Is3 = I3_Final(1:length(tspan)); % somatic input 
Id3 = I2_Final(1:length(tspan)); % dendritic input
%   paradigm=4; suprathreshold dendritic stimulation only
load('FiguresData\BACfiring\voltageIh4.mat', 'Vs', 'Vd');
load('FiguresData\BACfiring\stimulusIh4.mat', 'I3_Final', 'I2_Final');
Vs4 = Vs(1:length(tspan))+65; % voltage at the soma
Vd4 = Vd(1:length(tspan))+55; % voltage at the dendrites
Is4 = I3_Final(1:length(tspan)); % somatic input 
Id4 = I2_Final(1:length(tspan)); % dendritic input

% Note: there is a voltage differnce of 10 mV between the dendritic and
% somatic resting potential. We add the resting potential to the simulated
% voltages for visualization porposes. 

%% -- loading PC image and experimental results figure
[X,cmap] = imread('PC_ExpFig\PC.tif');
% experimental results (Schaefer et al., 2003)
[X2,cmap2] = imread('PC_ExpFig\Schaefer 2003b.tif');

%% -- creating the figure
figure('units','inch','position',[0,0,6.5,5]); % defining the figure's dimensions
scaleS = 25; % scale for the dendritic input
scaleD = 15; % scale for the somatic input
thre = 5.241203e-07; % [mA] input current threshold for Ca2+-spike
font = 10; % font size

%============== A
pos = [0.05 0 0.3 1]; % position of the plot
subplot('Position',pos)
imshow(X,cmap)
text(-150, 140, 'A', 'HorizontalAlignment','center', 'fontsize', font+4)
%============== B
pos = [0.38 0.81 0.3 0.25]; % position of the plot
subplot('Position',pos)
plot(tspan, Vs1, 'k', 'LineWidth', 1); % somatic voltage
hold on; 
plot(tspan, Vd1, 'r', 'LineWidth', 1, 'LineWidth', 1); % dendritic voltage
hold on
plot(tspan,10.*Is1./max(abs(Is1))-30,'k', 'LineWidth', 1); % somatic input current
hold on
plot(tspan,10.*Id1./max(abs(Id1))-30,'r', 'LineWidth', 1); % dendritic input current
hold on
plot(tspan,10.*ones(size(Id1)).*(thre)./max(abs(Id1))-30,'--k', 'LineWidth', 1); % dendritic threshold
hold off
text(0, 45, 'B', 'HorizontalAlignment','center', 'fontsize', font+4) % figure letter
text(mean(tspan), 45, 'Model', 'HorizontalAlignment','center', 'fontsize', font+2,'fontweight','bold') 
text(0, min(Vs1)+11, 'V_m', 'HorizontalAlignment','left', 'fontsize', font-1)
text(0, min(10.*Id1./max(abs(Id1))-8), 'I_{stim}', 'HorizontalAlignment','left', 'fontsize', font-1)
ylim([-30 80])
set(gca, 'Visible', 'off')
%============== C
pos = [0.38 0.53 0.3 0.25]; % position of the plot
subplot('Position',pos)
plot(tspan, Vs2, 'k', 'LineWidth', 1);  % somatic voltage
hold on; 
plot(tspan, Vd2, 'r', 'LineWidth', 1); % dendritic voltage
hold on
plot(tspan,20.*Is2./max(abs(Is2))-45,'k', 'LineWidth', 1); % somatic input current
plot(tspan,scaleD.*Id2./max(abs(Id2))-45,'r', 'LineWidth', 1); % dendritic input current
plot([60; 60], [40; 100], '-k',...
    [60; 70], [-7; -7], '-k', 'LineWidth', 1.5) % voltage and time scale
hold on
plot(tspan,scaleD.*ones(size(Id1)).*(thre)./max(abs(Id2))-45,'--k', 'LineWidth', 1); 
hold off
text(63, 80, '60 mV', 'HorizontalAlignment','left', 'fontsize', font+1) % voltage scale label
text(63, 50, '3 nA', 'HorizontalAlignment','left', 'fontsize', font+1) % current scale label
text(65, -20, '10 ms', 'HorizontalAlignment','center', 'fontsize', font+1) % time scale label
text(0, 100, 'C', 'HorizontalAlignment','center', 'fontsize', font+4) % figure letter
text(0, min(Vd2)+20, 'V_m', 'HorizontalAlignment','left', 'fontsize', font-1)
text(0, min(scaleS.*Is2./max(abs(Is2))-45)+15, 'I_{stim}', 'HorizontalAlignment','left', 'fontsize', font-1)
set(gca, 'Visible', 'off')
%============== D
pos = [0.38 0.27 0.3 0.25]; % position of the plot
subplot('Position',pos)
plot(tspan, Vs3, 'k', 'LineWidth', 1);  % somatic voltage
hold on; 
plot(tspan, Vd3, 'r', 'LineWidth', 1); % dendritic voltage
hold on
plot(tspan,scaleS.*Is3./max(abs(Is3))-50,'k', 'LineWidth', 1); % somatic input current
plot(tspan,scaleD.*Id3./max(abs(Id3))-90,'r', 'LineWidth', 1); % dendritic input current
hold on
plot(tspan,scaleD.*ones(size(Id1)).*(thre)./max(abs(Id3))-90,'--k', 'LineWidth', 1); % dendritic threshold 
text(0, 94, 'D', 'HorizontalAlignment','center', 'fontsize', font+4)
text(0, min(Vd3)+50, 'V_m', 'HorizontalAlignment','left', 'fontsize', font-1)
text(0, min(scaleS.*Is3./max(abs(Is3))-50)+23, 'I_{stim}', 'HorizontalAlignment','left', 'fontsize', font-1)
set(gca, 'Visible', 'off')
%============== E
pos = [0.38 0.02 0.3 0.25]; % position of the plot
subplot('Position',pos)
plot(tspan, Vs4, 'k', 'LineWidth', 1);  % somatic voltage
hold on; 
plot(tspan, Vd4, 'r', 'LineWidth', 1); % dendritic voltage
hold on
plot(tspan,scaleS.*Is4./max(abs(Is4))-45,'k', 'LineWidth', 1); % somatic input current
plot(tspan,scaleD.*Id4./max(abs(Id4))-45,'r', 'LineWidth', 1); % dendritic input current
hold on
plot(tspan,scaleD.*ones(size(Id1)).*(thre)./max(abs(Id4))-45,'--k', 'LineWidth', 1);  % dendritic threshold
text(0, 94, 'E', 'HorizontalAlignment','center', 'fontsize', font+4) % figure letter
text(0, min(Vd4)+45, 'V_m', 'HorizontalAlignment','left', 'fontsize', font-1)
text(0, min(scaleS.*Id4./max(abs(Id4))-50)+39, 'I_{stim}', 'HorizontalAlignment','left', 'fontsize', font-1)
set(gca, 'Visible', 'off')
legend({'Soma','Tuft'},'Location','none','position',[0.6 0.2 0.05 0.05], 'fontsize', font)
legend('boxoff')

% --- Experiment
pos = [0.7 0.005 0.3 1]; % position of the plot
subplot('Position',pos)
imshow(X2,cmap2)
annotation('textbox',[.77 .81 .1 .2],'String','Experiment','EdgeColor','none', 'fontsize', font+2,'fontweight','bold')

%% -- saving the figure
file = 'Figures';
if ~exist(file, 'dir') % checks if the folder already exists
    mkdir(file);  % creates a folder named 'file'
end

print(gcf, '-dtiff', '-r1000', 'Figures\fig_BACfiringExp.tiff');
