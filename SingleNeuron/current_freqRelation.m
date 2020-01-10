function current_freqRelation(V, step, tf, tspan, mustart_s, dmu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% current_freqRelation: computes the f-I relationship 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu = (mustart_s:dmu:(mustart_s+(round(tf/step))*dmu)); % [mA] mean stimulation current
mu = mu.*1e6; % convert mA-nA

intervals = 0:step:tf; % time intervals when the stimulation current had constant mean

mu = mu(1:length(intervals));

% detecting the APs generated at the soma during the simulation
[~, locs] = findpeaks(V, 'MinPeakProminence', 40, 'MinPeakDistance',6);
loc = tspan(locs); % times when the APs occurred 

freq = zeros(size(mu)); % pre-allocating memory, AP frequency per current step

for ii=1:(numel(freq)-1)
    ind = find(loc>=intervals(ii) & loc<intervals(ii+1)); % selecting the APs that occorred at each current step
    freq(ii) = numel(ind)./step; % computing the AP rate per current step
end

% fitting the f-I curve obtained with a linear function. Only frequencies
% different from zero were considered
ind = ismember(freq, 0); 
freq1 = freq(~ind); 
mu1 = mu(~ind);

P = polyfit(mu1,freq1,1);
yfit = P(1)*mu1 + P(2);

%% visualization
figure;
plot(mu(1:length(mu)-1).*1e3, freq(1:length(mu)-1), '.r', 'markersize', 18);
hold on;
plot(mu1, yfit, '-r');
ylabel('Mean spike rate (Ap/s)');
xlabel('Mean current (pA)');
set(gca,'linewidth',1,'fontsize',14,'fontweight','bold')
title(sprintf('Slope: %d',P(1)),'fontsize',18);

end