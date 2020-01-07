%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% fig_fIcurve
% Authors: B. Herrera et al. 2020
% This program generates the Figure 2. Frequency - Input relationship
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
clc

%% -- simulation time
t0 = 0; %[s] start time
tf = 24000*1e-3;  %[s] end time 
dt = 0.000001;  %[s] time increment

tspan = t0:dt:tf;   %[s] time span
tspan = tspan(1:(numel(tspan)-1));

%% Dendritic stimulation
% == loading simulated data
load('FiguresData\CurrentRamp\voltageIhRampD.mat', 'VsR'); % somatic voltage -> [time pts x trials]
load('FiguresData\CurrentRamp\stimulusIhRampD.mat', 'I2_Final', 'mean_u');
VsIhD = VsR;
VsIhPlotD = VsIhD(:,5)'; % selecting one trial to plot
mean_u2 = mean_u; % mean input current
clear mean_u VdR
clear VsR 

% == computing the f-I curve
stepD = 2; % [s] dendritic time step (total time with constant mu)
mustart_s = 0.2e-6; % [mA] mu start
dmu = 0.05e-6; % [mA] delta mu
muD = (mustart_s:dmu:(mustart_s+(round(tf/stepD))*dmu)); % [mA] dendritic mean current input
muD = muD.*1e6; % convert from mA->nA

% pre-allocating
PIhtD = zeros(2,length(VsIhD(1,:))); % linear fit

% time intervals with constant mean current
intervalsD = 0:stepD:tf; 
% selecting the mu values used for dendritic stimulation
muD = muD(1:length(intervalsD));

% pre-allocating frequency vector
freqIhD = zeros(length(muD)-1,length(VsIhD(1,:)));

for ii=1:length(VsIhD(1,:))
    % --- Detecting the spikes in each current step
    [~, locsIhD] = findpeaks(VsIhD(:,ii), 'MinPeakProminence', 40, 'MinPeakDistance',6);
    
    if ~isempty(locsIhD)
        locIhD = tspan(locsIhD);
    end
    clear locsIhD
    
    % --- Counting the number of spikes in each current step
    for jj=1:(numel(muD)-1)
        
        clear indIhD
       
        indIhD = find(locIhD>=intervalsD(jj) & locIhD<intervalsD(jj+1));
        freqIhD(jj,ii) = numel(indIhD)./stepD;
    end
    clear locIhD
    
    % --- linear fit of each trial
    freqtIhD = freqIhD(freqIhD(:,ii)>0,ii)'; %(5:end-1)
    mutIhD = muD(freqIhD(:,ii)>0);
    PIhtD(:,ii) = polyfit(mutIhD,freqtIhD,1);
end

% -- standard error of the mean frequency
stdE_PIhtD = std(PIhtD(1,:),0,2)./sqrt(length(PIhtD(1,:)));

clear VsIhD

%% Somatic stimulation
% == loading the simulated data
load('FiguresData\CurrentRamp\voltageIhRamp.mat', 'VsR');
load('FiguresData\CurrentRamp\stimulusIhRamp.mat', 'I3_Final', 'mean_u');
VsIh = VsR;
VsIhPlot = VsIh(:,5)'; % selecting one trial to plot
clear VsR

% == computing the f-I curve
step = 2; % [s] somatic time step (total time with constant mu)
mu = (mustart_s:dmu:(mustart_s+(round(tf/step))*dmu)); % [mA] somatic mean current input
mu = mu.*1e6; % convert from mA->nA

% pre-allocating
PIht = zeros(2,length(VsIh(1,:))); % linear fit

% time intervals with constant mean current
intervals = 0:step:tf; % somatic stimulation

% selecting the mu values used for somatic stimulation
mu = mu(1:length(intervals));

% pre-allocating frequency vector
freqIh = zeros(length(mu)-1,length(VsIh(1,:)));

for ii=1:length(VsIh(1,:))
    % --- Detecting the spikes in each current step
    [~, locsIh] = findpeaks(VsIh(:,ii), 'MinPeakProminence', 40, 'MinPeakDistance',6);
    
    if ~isempty(locsIh)
        locIh = tspan(locsIh);
    end
    clear locsIh
    
    % --- Counting the number of spikes in each current step
    for jj=1:(numel(mu)-1)
        
        clear indIh
        
        indIh = find(locIh>=intervals(jj) & locIh<intervals(jj+1));
        freqIh(jj,ii) = numel(indIh)./step;
    end
    clear locIh
    
    % --- linear fit of each trial
    freqtIh = freqIh(freqIh(:,ii)>0,ii)'; %(5:end-1)
    mutIh = mu(freqIh(:,ii)>0);
    PIht(:,ii) = polyfit(mutIh,freqtIh,1);
end

% --- standard error of the mean frequency
stdE_PIht = std(PIht(1,:),0,2)./sqrt(length(PIht(1,:)));

clear VsIh

%% --- mean rate predicted by the model for somatic and dendritic stimulations, respectively
% - Somatic stimulation
freqIhM = mean(freqIh,2); % mean frequency
stdE_freqIhKsM = std(freqIh,0,2)./sqrt(length(freqIh(1,:))); % standard error of the mean frequency
% - Dendritic stimulation
freqIhMD = mean(freqIhD,2); % mean frequency
stdE_freqIhKsMD = std(freqIhD,0,2)./sqrt(length(freqIhD(1,:))); % standard error of the mean frequency

%% --- Computing the slope of the f-I curves
% - Somatic stimulation
freq1Ih = freqIhM';
mu1Ih = mu(1:end-1);
PIh = polyfit(mu1Ih,freq1Ih,1); % Linear fit
yfitIh = PIh(1)*mu1Ih + PIh(2);
% R^2, Goodness of Fit 
RsqFitIh = 1 - sum((freq1Ih - yfitIh).^2)/sum((freq1Ih - mean(freq1Ih)).^2); 
% - Dendritic stimulation
freq1IhD = freqIhMD(freqIhMD>1)';
mu1IhD = muD(freqIhMD>1);
PIhD = polyfit(mu1IhD,freq1IhD,1);
yfitIhD = PIhD(1)*mu1IhD + PIhD(2);
% R^2, Goodness of Fit 
RsqFitIhD = 1 - sum((freq1IhD - yfitIhD).^2)/sum((freq1IhD - mean(freq1IhD)).^2); 

%% --- Experimental f-I curves
% -- Somatic stimulation
% - Larkum et al. 2004 results
freqLarkum = [1, 2, 3.9, 6, 7, 8, 12.5, 14.5, 17, 20, 22, 24]; % mean AP rate
ILarkum = 1e-3.*(200:50:750); % [nA] mean input current
PL = polyfit(ILarkum,freqLarkum,1); % Linear Fit
yfitL = PL(1)*ILarkum + PL(2);
RsqFitL = 1 - sum((freqLarkum - yfitL).^2)/sum((freqLarkum - mean(freqLarkum)).^2); % R^2, Goodness of Fit 
% - Bahl et al. 2012 results
freqBahl = [2 8 10 14 16 18 19 21 22 26 27 28]; % 30 31 32 36 38];  % mean AP rate
IBahl = 1e-3.*(200:50:750); % [nA] mean input current
PBahl = polyfit(IBahl,freqBahl,1);  % Linear Fit
yfitBahl = PBahl(1)*IBahl + PBahl(2);
% - Mean slope between both experimental results
yfitExpMean = mean([yfitL; yfitBahl],1);

% -- Dendritic stimulation
% - Larkum et al. 2004 results
freqLarkumD = [0 2 4 7 8 10 12 14 16]; % mean AP rate
ILarkumD = 1e-3.*(550:50:950); % [nA] mean input current

%% --- Distance between somatic and dendritic stimulation's curves
% - Larkum et al. 2004 results
deltaILarkum = ILarkumD(ismember(freqLarkumD, [2 4 7 8 12 14]))- ...
    ILarkum(ismember(freqLarkum, [2 3.9 7 8 12.5 14.5])); % difference
deltaILarkumMean = mean(deltaILarkum); % mean difference
deltaILarkumStd = std(deltaILarkum); % standard deviation
% - Simulated data
freq = freq1Ih(1:6);
predmuD = (freq-PIhD(2))./PIhD(1);
deltaImodel = predmuD-mu1Ih([1 2 3 4 5 6]); % difference
deltaImodelMean = mean(deltaImodel); % mean difference
deltaImodelStd = std(deltaImodel); % standard deviation

% - t-test between the observed and predicted distances
[h, p,CI,STATS] = ttest(deltaILarkum, deltaImodel);

%% -- loading the PC image
[X,cmap] = imread('PC_ExpFig\PC_den.tif');

%% -- creating the figure
figure('units','inch','position',[0,0,6,6]); % defining the dimensions of the figure
font = 10; % font size

tfplot = 16; % [s] time span for voltage plots

% ================ A
% --- PC image
pos = [0 0.55 0.3 0.45]; % position of the plot
subplot('Position',pos)
imshow(X,cmap)
text(-100, 80, 'A', 'HorizontalAlignment','center', 'fontsize', font+4,'fontweight','bold')

% ================ B
% --- Somatic voltage for somatic stimulation
pos1 = [0.3 0.88 0.685 0.115]; % position of the plot [0.3 0.64 0.685 0.125]
subplot('Position',pos1)
plot(tspan((0*round(1/dt)+1):(tfplot*round(1/dt)+1)), VsIhPlot((0*round(1/dt)+1):(tfplot*round(1/dt)+1)),...
    '-b', 'LineWidth', 0.5); 
ylim([min(VsIhPlot) max(VsIhPlot)])
set(gca, 'Visible', 'off')
text(-0.05, 30, 'B', 'HorizontalAlignment','right', 'fontsize', font+4,'fontweight','bold') % letter of the figure

% --- Somatic input current
pos2 = [0.3 0.78 0.685 0.1]; % position of the plot
subplot('Position',pos2)
plot(tspan((0*round(1/dt)+1):(tfplot*round(1/dt)+1)), 1e6.*I3_Final((0*round(1/dt)+1):(tfplot*round(1/dt)+1)),...
    '-b', 'LineWidth', 1); % plotting the current
hold on;
plot(tspan((0*round(1/dt)+1):(tfplot*round(1/dt)+1)), 1e6.*mean_u((0*round(1/dt)+1):(tfplot*round(1/dt)+1)),...
    '-k', 'LineWidth', 0.5); hold on
set(gca, 'Visible', 'off')

% --- Somatic voltage for dendritic stimulation
pos1 = [0.3 0.66 0.685 0.115]; % position of the plot
subplot('Position',pos1)
plot(tspan((0*round(1/dt)+1):(tfplot*round(1/dt)+1)), VsIhPlotD((0*round(1/dt)+1):(tfplot*round(1/dt)+1)),...
    '-b', 'LineWidth', 0.5); 
hold on
plot([4; 4], [-30; 20], '-k', 'LineWidth', 1.5) % voltage scale
hold on
plot([tspan((1*round(1/dt)+1)); tspan((3*round(1/dt)+1))], [-80; -80], '-k', 'LineWidth', 1.5) % time scale
hold off
text(3.8, mean([-30 20]), '50 mV', 'HorizontalAlignment','right', 'fontsize', font) % voltage scale label
text(mean([tspan((1*round(1/dt)+1)) tspan((3*round(1/dt)+1))]), -100, '2 s', 'HorizontalAlignment','center', 'fontsize', font) % time scale label
ylim([min(VsIhPlot) max(VsIhPlot)])
set(gca, 'Visible', 'off')

% --- Dendritic input current
pos2 = [0.3 0.54 0.685 0.1]; % position of the plot
subplot('Position',pos2)
plot(tspan((0*round(1/dt)+1):(tfplot*round(1/dt)+1)), 1e6.*I2_Final((0*round(1/dt)+1):(tfplot*round(1/dt)+1)), '-r', 'LineWidth', 1);
hold on;
plot(tspan((0*round(1/dt)+1):(tfplot*round(1/dt)+1)), 1e6.*mean_u2((0*round(1/dt)+1):(tfplot*round(1/dt)+1)),'-k', 'LineWidth', 0.5); hold on;
plot([10; 12], [-0.1; -0.1], '-k', 'LineWidth', 1.5);  % time scale
hold on
plot([12; 12], [-0.1; 0.1], '-k', 'LineWidth', 1.5) % current scale
hold off
text(mean([10 12]), -0.3, '2 s', 'HorizontalAlignment','center', 'fontsize', font) % time scale label
text(12.1, 0, '0.2 nA', 'HorizontalAlignment','left', 'fontsize', font) % current scale label
set(gca, 'Visible', 'off')

% ================ C
% --- f-I curve
pos3 = [0.1 0.08 0.45 0.45]; % position of the plot
subplot('Position',pos3)
p1 = plot(ILarkum, freqLarkum, 'k', 'linewidth',1.5); % experimental curve by Larkum et al. 2004
hold on; 
plot(IBahl, freqBahl, 'k', 'linewidth',1.5); % experimental curve by Bahl et al. 2012
% filling the area between these curve
x2 = [ILarkum, fliplr(ILarkum)];
inBetween = [freqLarkum, fliplr(freqBahl)];
fill(x2, inBetween, 1.8.*[0.5,0.5,0.5]);
hold on
p2 = plot(ILarkum, yfitExpMean, '--k', 'linewidth',1.5); % mean linear fit - experiments
hold on; 
p3 = errorbar(mu1Ih, freq1Ih, 1.96.*stdE_freqIhKsM, '.b', 'markersize', 8); % simulated f-I for somatic stimulation
hold on
p4 = plot(mu1Ih, yfitIh, '--b', 'linewidth',1.5); % linear fit of the f-I for somatic stimulation
hold on; 
p5 = errorbar(muD(1:end-1), freqIhMD, 1.96.*stdE_freqIhKsMD, '.r', 'markersize', 8);  % simulated f-I for dendritic stimulation
hold on; 
p6 = plot(mu1IhD, yfitIhD, '--r', 'linewidth',1.5); % linear fit of the f-I for dendritic stimulation
hold on
plot([mu1Ih(4), mu1IhD(6)], [7.12 7.12], '-k', 'LineWidth',1) % delta I 
text(mean([mu1Ih(4), mu1IhD(6)]), 6, '\DeltaI', 'HorizontalAlignment','center', 'fontsize', font)
hold off
ylabel('Mean Spike Rate (AP/s)');
xlabel('Mean Current (nA)');
legend([p1 p2 p3 p4 p5 p6], 'Experiment', sprintf('Fit (Slope: %4.2f, R^2: %4.3f)',PL(1),round(RsqFitL,3)),...
    'Model - Somatic Input',...
    sprintf('Fit (Slope: %4.2f, R^2: %4.3f)',PIh(1),round(RsqFitIh,3)),...
    'Model - Distal Apical Input ',...
    sprintf('Fit (Slope: %4.2f, R^2: %4.3f)',PIhD(1),round(RsqFitIhD,3)),...
    'Location','northoutside', 'fontsize', 6,'NumColumns',2)
legend('boxoff')
set(gca, 'box', 'off','linewidth',1,'fontsize',font,'fontweight','bold')
xlim([0.2 0.76])
% letter of the figure
annotation('textbox','position',[0 0.31 .25 .25],'String','C','fontsize', font+4,'fontweight','bold','EdgeColor','none')

% ================ D
pos3 = [0.65 0.15 0.3 0.25]; % position of the plot
subplot('Position',pos3)
Xb = [1 3]; % position of the bars
Y = [deltaILarkumMean, deltaImodelMean]; % mean delta I values
p1 = bar(Xb(1),Y(1));
hold on;
p2 = bar(Xb(2),Y(2));
set(p1,'FaceColor',[0.5,0.5,0.5]);
set(p2,'FaceColor','blue');
hold on
errorbar(Xb(1), Y(1), deltaILarkumStd, 'Color', 'k')
hold on
errorbar(Xb(2), Y(2), deltaImodelStd, 'Color', 'k')
hold on
plot([Xb(1)-0.015, Xb(2)+0.015], [1 1].*max(Y(1,1:2))+0.1, '-k', 'LineWidth',2)
if p<0.05 % if the difference between the observed and simulated deltaI is significant
    plot(mean([Xb(1)-0.015, Xb(2)+0.015]), max(Y(1,1:2))+0.12, '*k')
else % if the difference between the observed and simulated deltaI is not significant
    text(mean([Xb(1)-0.015, Xb(2)+0.015]), max(Y(1,1:2))+0.13, 'NS', 'HorizontalAlignment','center', 'fontsize', font-1)
end
hold off
ylim([0 0.5])
xticks(Xb)
xticklabels({'Experiment', 'Model'})
ylabel('\DeltaI (nA)')
set(gca, 'box', 'off','linewidth',1,'fontsize',font,'fontweight','bold')
text(-1, 0.6, 'D', 'HorizontalAlignment','right', 'fontsize', font+4,'fontweight','bold') % letter of the figure

%% -- saving the figure
file = 'Figures';
if ~exist(file, 'dir') % checks if the folder already exists
    mkdir(file);  % creates a folder named 'file'
end

print(gcf, '-dtiff', '-r1000', 'Figures\fig_fIcurve.tiff');
