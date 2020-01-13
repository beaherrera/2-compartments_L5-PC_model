%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% fig_CSD
% Authors: B. Herrera et al. 2020
% This program generates the Figure 5. LFP reflections of the Dendritic 
% Ca2+-spikes – CSD analysis. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%% simulation time
t0 = 0; %[s] start time
tf = 80*1e-3;  %[s] end time
dt = 0.000001;  %[s] time increment

ts = t0:dt:tf;   %[s] time span
ts = ts(1:(numel(ts)-1));

% somatic stimulation window
I3_onset = 10e-3; % stimulus onset
I3_end = 30e-3; % stimulus end

%% -- loading data

% LFP
file = 'FiguresData\data1000n';
load([file '\lfp_L5PCs.mat'], 'lfp');
Fs = 10e3; % [Hz] sampling frequency
Ve = eegfilt(lfp, Fs, 0, 90,0,250); % filtering the LFP at 90 Hz

% CSD
% no Ih
load([file '\CSDNoIh.mat'], 'zs', 'iCSD');
zsNoIh = zs;
iCSDnoIh = iCSD;
clear zs iCSD
% Ih
load([file '\CSD.mat'], 'tspan', 'zs', 'iCSD');

% ICaL
load('FiguresData\data1000n\ICaLTrial1.mat'); 
ICaLIh = mean(ICaL,2);
clear ICaL
load('FiguresData\data1000n\ICaLNoIhTrial1.mat'); 
ICaLnoIh =  mean(ICaL,2);
clear ICaL

%% detecting spikes for raster plot and the histogram

trial=1; % selecting a trial

load([file '\voltageTrial' num2str(trial) '.mat'], 'voltage'); % loading the L5-PCs' voltages of the selected trial
Vsoma = voltage.Vs; % somatic voltage for all neurons
Vd = voltage.Vd; % dendritic voltage for all neurons
clear voltage

for ii=1:length(Vsoma(1,:)) % loop for the number of neurons
    [~, locsSoma] = findpeaks(Vsoma(:,ii), 'MinPeakProminence', 50, 'MinPeakDistance',6); % detecting the somatic APs
    [~, locsD] = findpeaks(Vd(:,ii), 'MinPeakProminence', 50, 'MinPeakDistance',6); % detecting the dendritic Ca spikes
    
    timesSoma = ts(locsSoma); % creates a vector with the spikes' time -> Soma
    timesD = ts(locsD); % creates a vector with the spikes' time -> Dendrites
    
    % creating a logic matrix of the spikes times for the raster plot
    if ii==1
        spk_timesSoma = ismember(ts, timesSoma);
        spk_timesD = ismember(ts, timesD);
    else
        spk_timesSoma = cat(1, spk_timesSoma, ismember(ts, timesSoma));
        spk_timesD = cat(1, spk_timesD, ismember(ts, timesD));
    end
    
    clear locsSoma locsD timesD timesSoma
end

% creating a vector with all the spikes' time of all neurons to construct
% the histogram
tsM = ts.*ones(length(Vsoma(1,:)),1);
Spks_S = tsM(spk_timesSoma); % Somatic spikes' time
Spks_D = tsM(spk_timesD); % Ca spikes' time

%% Ih-noIh comparison

% electrodes' location
Ne = 16; % number of electrodes in the shank
a = 0.1; % [mm] position of the first electrode
elec_spacing = 0.1; % [mm] electrode spacing
ze = a:elec_spacing:((Ne-1)*elec_spacing + a); % [mm] electrode positions with respect to the pia surface

% loading the CSD for trials with Ih
load([file '\CSDT.mat'], 'zsT', 'iCSDT');
zsTIh = zsT;
iCSDTIh = iCSDT;
clear zsT iCSDT
% loading the CSD for trials without Ih
load([file '\CSDNoIhT.mat'], 'zsT', 'iCSDT');
zsTnoIh = zsT;
iCSDTnoIh = iCSDT;
clear zsT iCSDT

% pre-allocating memory
% - Ih present
CSDIh_top = zeros(10,length(iCSDTIh(1,:,1))); % CSD at 0.4 mm below the pia matter
iAmpIh = zeros(10,1); % index of the maximum value of the CSD at 0.4 mm below the pia matter
maxAmpIh = zeros(10,1); % maximum value of the CSD at 0.4 mm below the pia matter
areaIh = zeros(10,1); % area of the CSD^2 trace at 0.4 mm below the pia matter
% - Ih blocked
CSDnoIh_top = zeros(10,length(iCSDTIh(1,:,1))); % CSD at 0.4 mm below the pia matter 
iAmpnoIh = zeros(10,1); % index of the maximum value of the CSD at 0.4 mm below the pia matter
maxAmpnoIh = zeros(10,1); % maximum value of the CSD at 0.4 mm below the pia matter 
areanoIh = zeros(10,1); % area of the CSD^2 trace at 0.4 mm below the pia matter 


for ii=1:length(iCSDTIh(1,1,:)) % loop for the trials
    % finding the index in the CSD that correspond to the depth: 0.4 mm 
    [~, indnoIhT] = min(abs(zsTnoIh(:,ii)-ze(4))); % Ih blocked
    [~, indIhT] = min(abs(zsTIh(:,ii)-ze(4))); % Ih present
    % finding the indexes of the maximum value of the CSD at 0.4 mm below the pia matter
    [~, iAmpnoIh(ii)] = max(abs(iCSDTnoIh(indnoIhT,:,ii)),[],2); % Ih blocked
    [~, iAmpIh(ii)] = max(abs(iCSDTIh(indIhT,:,ii)),[],2); % Ih present
    
    % Ih blocked
    maxAmpnoIh(ii,1) = iCSDTnoIh(indnoIhT,iAmpnoIh(ii),ii); % maximum value of the CSD at each depth
    CSDnoIh_top(ii,:) = iCSDTnoIh(indnoIhT,:,ii); % CSD at each depth
    areanoIh(ii,1) = trapz(tspan, squeeze(CSDnoIh_top(ii,:)).^2); % area of the CSD^2 trace at each depth
    % Ih present
    maxAmpIh(ii,1) = iCSDTIh(indIhT,iAmpIh(ii),ii); % maximum value of the CSD at each depth
    CSDIh_top(ii,:) = iCSDTIh(indIhT,:,ii); % CSD at each depth
    areaIh(ii,1) = trapz(tspan, squeeze(CSDIh_top(ii,:)).^2); % area of the CSD^2 trace at each depth
end

% Wilcoxon signed rank test to compare the differnce between the sink at
% 0.4 below the pia matter with and without Ih. The comparison was done
% using the area below the square of the CSD traces at that depth 
[p,h] = signrank(areaIh,areanoIh);

%% visualization

figure('units','inch','position',[0,0,6.5,7]); % defining the size of the figure
font = 10; % font size

% ====================== A
% =================== one neuron response
Vs1n = Vsoma(:,20);
Vd1n = Vd(:,20);

% --- Somatic voltage
pos1 = [0.025 0.82 0.25 0.15]; % position of the plot
subplot('Position',pos1)
% plotting the voltage
plot(ts, Vs1n, '-k', 'LineWidth', 1);
hold on
plot(ts(end)-[0; 10e-3], [-80; -80], '-k', 'LineWidth', 1.5);  % time scale
hold on
plot(ts(end)-[10e-3; 10e-3], [-50; 0], '-k', 'LineWidth', 1.5);  % voltage scale
hold off
text(ts(end)-mean([0 10e-3]), -95, '10 ms', 'HorizontalAlignment','center', 'fontsize', font) % time scale label
text(ts(end)-8e-3, mean([-50; 0]), '50 mV', 'HorizontalAlignment','left', 'fontsize', font) % voltage scale label
text(5e-3, -30, 'V_s', 'HorizontalAlignment','right', 'fontsize', font,'Color','k')
text(0, 60, 'A', 'HorizontalAlignment','right', 'fontsize', font+4,'fontweight','bold') % letter of the figure
set(gca, 'Visible', 'off')

% --- dendritic voltage
pos2 = [0.025 0.65 0.25 0.15]; % position of the plot
subplot('Position',pos2)
plot(ts, Vd1n, '-r', 'LineWidth', 1);
hold on
plot([I3_onset; I3_end], [-75; -75], '-b', 'LineWidth', 4.5);  % stimulation window
text(5e-3, -30, 'V_d', 'HorizontalAlignment','right', 'fontsize', font,'Color','red')
set(gca, 'Visible', 'off')

% ====================== B
% ============================ raster plot and histogram
pos1 = [0.39 0.78 0.26 0.2]; % position of the plot
subplot('Position',pos1)
% raster plot of the APs fired by 100 L5-PCs in the network 
LineFormat = struct();
LineFormat.Color = [0 0 0]; % color of the vertical lines
LineFormat.LineWidth = 1.5; % width of the vertical lines
plotSpikeRaster(spk_timesSoma(1:100,:),'PlotType','vertline','LineFormat',LineFormat,'VertSpikeHeight',5);
hold on
% raster plot of the Ca2+ spikes fired by 100 L5-PCs in the network 
LineFormat = struct();
LineFormat.Color = [1 0 0]; % color of the vertical lines
LineFormat.LineWidth = 1.5; % width of the vertical lines
plotSpikeRaster(spk_timesD(1:100,:),'PlotType','vertline','LineFormat',LineFormat,'VertSpikeHeight',5);
hold off
dim2 = [.418 .785 .062 .002];
annotation('rectangle',dim2,'Color','b','FaceColor','b') % stimulation window
ylabel('L5-PCs')
xticks([0 2 4 6 8].*1e4)
xticklabels({})
ylim([0 105])
text(-1e4, 0, 'B', 'HorizontalAlignment','right', 'fontsize', font+4,'fontweight','bold') % letter of the figure
box 'off'
set(gca,'linewidth',1.5,'fontsize',font-1,'fontweight','bold')

% --- Histogram
pos1 = [0.39 0.66 0.26 0.1]; % position of the plot
subplot('Position',pos1)
hist1 = histogram(Spks_S.*1e3); % somatic APs
hold on
hist2 = histogram(Spks_D.*1e3,'FaceAlpha',1); % Ca2+ spikes
hist1.FaceColor = 'k';
hist2.FaceColor = 'r';
box 'off'
xlabel('Time (ms)')
ylabel('Events')
xticks([0 20 40 60 80])
set(gca,'linewidth',1.5,'fontsize',font-1,'fontweight','bold')
text(80, 250, 'Tuft', 'HorizontalAlignment','right', 'fontsize', font-1,'fontweight','bold', 'color', 'r')
text(80, 200, 'Soma', 'HorizontalAlignment','right', 'fontsize', font-1,'fontweight','bold', 'color', 'k')

% ====================== C
% ======================================== LFP
pos1 = [0.755 0.58 0.2 0.4]; % position of the plot
subplot('Position',pos1)
LFPmax = 7*max(abs(lfp),[],'all'); % scale for the LFP plot
plot(tspan, (ze.*ones(length(tspan), length(ze)))','-', 'Color', [0.5 0.5 0.5], 'linewidth', 1)
hold on
for ii = 1:Ne
    plot(tspan,Ve(Ne-ii+1,:)./LFPmax + ze(ii), ...
        'color', 'r', 'clipping','on', 'linewidth', 0.5)
end
hold on
plot(I3_onset*ones(length(ze)+2).*1e3, [0 ze 1.7], '--', 'Color', [0.5 0.5 0.5], 'linewidth', 1)
hold on
plot([I3_onset; I3_end].*1e3, [0; 0], '-b', 'LineWidth', 2); % stimulation window
hold on
plot(tspan(end)-[20 40], [0; 0], '-k', 'LineWidth', 1.5);  % time scale
hold on
plot([tspan(end) tspan(end)]-20, [-0.1 ; -0.05]+ ze(1), '-k', 'LineWidth', 1.5); % voltage scale
hold off
text(tspan(end)-mean([20 40]), -0.05, '20 ms', 'HorizontalAlignment','center', 'fontsize', font) % time scale label
text(tspan(end)-15, mean([-0.1 ; -0.05]+ ze(1)), '0.5 mV', 'HorizontalAlignment','left', 'fontsize', font) % voltage scale label
ylabel('Depth (mm)','fontsize',font,'fontweight','bold')
yticks(0.2:0.2:1.6)
tickslabLFP = {'1.5'};
for ii=1.3:-0.2:0.1
    tickslabLFP = cat(2, tickslabLFP, num2str(ii,'%4.1f'));
end
yticklabels(tickslabLFP)
ax = gca;
ax.FontSize = font;
ylim([0 1.7])
box 'off'
ax.YRuler.Axle.LineStyle = 'none';
set(ax,'XTickLabel',[],'XTick',[],'XColor',[1 1 1])
text(-15000e-3, 1.7, 'C', 'HorizontalAlignment','right', 'fontsize', font+4,'fontweight','bold') % letter of the figure

% ====================== D
% =========================================== CSD
max_CSD = max(abs([iCSDnoIh iCSD]), [], 'all');
iCSDmax = 10*max(abs(iCSD),[],'all'); % scale for the CSD traces plot

% ===================== noIh
pos1 = [0.08 0.03 0.2 0.5]; % position of the plot
subplot('Position',pos1)
ind1 = find(zs>=0.09 & zs<=0.1);
ind2 = find(zs>=1.59 & zs<=1.6);
imagesc(tspan, 1.7-zsNoIh(ind1:ind2), iCSDnoIh(ind1:ind2,1:length(tspan)));
hold on
% CSD traces
for ii = 1:Ne
    [~, indnoIh] = min(abs(zs-ze(ii)));
    plot(tspan, (iCSDnoIh(indnoIh,:)./iCSDmax +1.7- zs(indnoIh)), ...
        'color', [128 0 32]./255, 'clipping','on', 'linewidth', 0.5)
end
hold on
plot(I3_onset*ones(length(ze)).*1e3, 1.7-ze, '--k', 'linewidth', 1)
hold on
plot(tspan(end)-[0 20], 1.7-[1.61; 1.61], '-k', 'LineWidth', 1.5);  % time scale
hold on
plot([I3_onset; I3_end].*1e3, 1.7-[1.61; 1.61], '-b', 'LineWidth', 2); % stimulation window
hold off
ylim([0.05 1.7])
text(tspan(end)-mean([0 20]), 1.7-1.67, '20 ms', 'HorizontalAlignment','center', 'fontsize', font) % time scale label
ax = gca;
ax.FontSize = font;
ax.YDir = 'normal';
colormap(jet);
bar_min = -max_CSD;
bar_max = max_CSD;
caxis([bar_min bar_max]);
ylabel('Depth (mm)','fontsize',font,'fontweight','bold')
yticks(0.2:0.2:1.6)
tickslab = {'1.5'};
for ii=1.3:-0.2:0.1
    tickslab = cat(2, tickslab, num2str(ii,'%4.1f'));
end
yticklabels(tickslab)
set(ax,'XTickLabel',[],'XTick',[],'XColor',[1 1 1])
ax.YRuler.Axle.LineStyle = 'none';
text(mean([tspan(1) tspan(end)]), 1.75-0, 'Blocking Ih', 'HorizontalAlignment','center', 'fontsize', font+4,'fontweight','bold')
text(-23000e-3, 1.75-0, 'D', 'HorizontalAlignment','right', 'fontsize', font+4,'fontweight','bold') % letter of the figure

% ========================= Ih
pos3 = [0.3 0.03 0.2 0.5]; % position of the plot
subplot('Position',pos3)
ind1 = find(zs>=0.09 & zs<=0.1);
ind2 = find(zs>=1.59 & zs<=1.6);
imagesc(tspan, 1.7-zs(ind1:ind2), iCSD(ind1:ind2,1:length(tspan)));
hold on
% CSD traces
for ii = 1:Ne
    [~, indIh] = min(abs(zs-ze(ii)));
    plot(tspan, (iCSD(indIh,:)./iCSDmax + 1.7-zs(indIh)), ...
        'color', [128 0 32]./255, 'clipping','on', 'linewidth', 0.5)
end
hold on
plot(I3_onset*ones(length(ze)).*1e3, ze, '--k', 'linewidth', 1)
hold on
plot([I3_onset; I3_end].*1e3, 1.7-[1.61; 1.61], '-b', 'LineWidth', 2); % stimulation window
hold on
plot([tspan(end) tspan(end)], ([0 ; 5])./iCSDmax + ze(1), '-k', 'LineWidth', 1.8); % CSD scale
hold off
text(tspan(end)+5, mean([0 ; 0.053]+ ze(1)), '5 \muA/mm^3', 'HorizontalAlignment','left', 'fontsize', font) % CSD scale label
ylim([0.05 1.7])
c = colorbar;
colormap(jet);
ax = gca;
ax.FontSize = font;
ax.YDir = 'normal';
c.Label.String = 'Amplitude (\muA/mm^{3})';
c.Label.FontSize = font;
c.FontSize = font;
c.FontWeight = 'bold';
c.Position = [0.51 0.1 0.02 0.2];
c.TickLength = 0.05;
bar_min = -max_CSD;
bar_max = max_CSD;
caxis([bar_min bar_max]);
yticks(0.1:0.2:1.5)
yticklabels({})
xtick = [];
box 'off'
set(ax,'XTickLabel',[],'XTick',[],'XColor',[1 1 1],'fontweight','bold')
ax.YRuler.Axle.LineStyle = 'none';
text(mean([tspan(1) tspan(end)]), 1.75-0, 'Ih', 'HorizontalAlignment','center', 'fontsize', font+4,'fontweight','bold')

% ============================= selecting zoomed plot
annotation('rectangle',[0.36 0.34 0.07 0.11],'Color','k');
h2 = annotation('line',[0.43 0.66],[0.45 0.55]);
h2.LineStyle = '--';
h2.LineWidth = 1.5;
h2.Color = [0.5 0.5 0.5];
h3 = annotation('line',[0.43 0.66],[0.34 0.36]);
h3.LineStyle = '--';
h3.LineWidth = 1.5;
h3.Color = [0.5 0.5 0.5];

% ============================= zoomed plot
pos1 = [0.66 0.36 0.2 0.19]; % position of the plot
subplot('Position',pos1)
Llim = 15;
Ulim = 55;
iCSDn = iCSD((zs>=0.3) & (zs<=0.6),:);
tspanS = tspan((tspan>=Llim) & (tspan<=Ulim))-10;
plot(tspanS, (ze(3:6).*ones(length(tspan((tspan>=Llim) & (tspan<=Ulim))), length(ze(3:6))))','--', 'Color', [0.5 0.5 0.5], ...
    'linewidth', 1)
hold on
jj = 6;
for ii = 3:6
    [~, indnoIh] = min(abs(zs((zs>=0.3) & (zs<=0.6))-ze(ii)));
    plot(tspanS, (iCSDn(indnoIh,(tspan>=Llim) & (tspan<=Ulim))./iCSDmax + ze(jj)), ...
        'color', [128 0 32]./255, 'clipping','on', 'linewidth', 1)
    indT = find((iCSDn(indnoIh,(tspan>=Llim) & (tspan<=Ulim))./iCSDmax + ze(jj))<ze(jj)-0.0028,1);
    hold on
    plot(tspanS(indT),iCSDn(indnoIh, indT)./iCSDmax + ze(jj)-0.028, '^k', 'MarkerFaceColor', 'k', 'markersize', 8)
    
    jj = jj-1;
end
hold on
plot([45 45], [0 ; 5]./iCSDmax+ 0.4, '-k', 'LineWidth', 1.8); % CSD scale
hold off
text(46.5, mean([0 ; 5]./iCSDmax + 0.4), '5 \muA/mm^3', 'HorizontalAlignment','left', 'fontsize', font) % CSD scale label
ylabel('Depth (mm)')
xlabel('Time after light onset (ms)')
yticks(0.3:0.1:0.6)
xticks([20 40])
tickslab = {'0.6'};
for ii=0.5:-0.1:0.3
    tickslab = cat(2, tickslab, num2str(ii,'%4.1f'));
end
yticklabels(tickslab)
ylim([0.25 0.65])
xlim([min(tspanS) max(tspanS)])
box 'off'
set(gca,'linewidth',1.5,'fontsize',font-1,'fontweight','bold')

% ====================================== ICaL
pos2 = [0.7 0.225 0.2 0.05]; % position of the plot
subplot('Position',pos2)
plot(ts.*1e3, ICaLIh.*1e-6, '-', 'Color',[105 105 105]./256, 'LineWidth', 1) % Ih present
hold on
plot(ts.*1e3, ICaLnoIh.*1e-6, '-', 'Color',[0 206 209]./256, 'LineWidth', 1) % Ih blocked
ylim([-7 0].*1e-11)
xlabel('Time (ms)')
ylabel('I_{CaL} (A)')
box 'off'
set(gca,'linewidth',1.5,'fontsize',8,'fontweight','bold')

% ====================================== Ih Vs noIh
pos3 = [0.7 0.05 0.2 0.118]; % position of the plot
subplot('Position',pos3)
X = 0.4;
Y = [mean(maxAmpIh)', mean(maxAmpnoIh)'];
stdErrorIh = std(maxAmpIh)./sqrt(length(maxAmpIh));
stdErrornoIh = std(maxAmpnoIh)./sqrt(length(maxAmpnoIh));
b = bar([0.3 X], Y.*ones(2,2)); % bar plot
b(1).FaceColor = [105 105 105]./256; % color for Ih present
b(2).FaceColor = [0 206 209]./256; % color for Ih blocked
hold on
errorbar(X+0.015, Y(2), stdErrornoIh, 'Color', 'k')
hold on
errorbar(X-0.015, Y(1), stdErrorIh, 'Color', 'k')
hold on
plot([X-0.015, X+0.015], [1 1], '-k', 'LineWidth',2)
if p<0.05 % if there is a significant difference
    plot(mean([X-0.015, X+0.015]), 1.75, '*k')
else % if there is not a significant difference
    text(mean([X-0.015, X+0.015]), 1.75, 'NS', 'HorizontalAlignment','center', 'fontsize', font-1)
end
hold off
xlim([0.35 0.45])
ylim([-6 1.75])
xticks(0.4)
xlabel('Depth (mm)')
ylabel('Amplitude (\muA/mm^{3})')
% legend
text(0.45, 2.5, 'Ih', 'HorizontalAlignment','left', 'fontsize', 8,'fontweight','bold', 'color', [105 105 105]./256)
text(0.45, 1.5, 'Blocking Ih', 'HorizontalAlignment','left', 'fontsize', 8,'fontweight','bold', 'color', [0 206 209]./256)
box 'off'
set(gca,'linewidth',1.5,'fontsize',8,'fontweight','bold')

%% -- saving the figure
file = 'Figures';
if ~exist(file, 'dir') % checks if the folder already exists
    mkdir(file);  % creates a folder named 'file'
end

print(gcf, '-dtiff', '-r1000', 'Figures\fig_CSD.tiff');

