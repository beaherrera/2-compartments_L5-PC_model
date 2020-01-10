function [I, mean_u] = Inoise_incMu(tspan, tf, dt, mustart_s, dmu, step, sigma)
%% Creates an in vivo-like input using the Ornstein-Uhlenbeck method. 
% Generates noisy currents with mean mu, standard deviation sigma, and
% correction length tau according to the equation:
% I(t+dt)=I(t)+(mu-I(t))/tau+sigma*Gt*sqrt(2*dt/tau)
% where: Gt is a random number taken each time step from a Gaussian
% distribution with mean and standard deviation 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

len = size(tspan); % number of time pts in the integration

% - current parameters
tau = 3e-3; % [s] time correlation length

% - current mean
mu = (mustart_s:dmu:(mustart_s+(round(tf/step))*dmu)); % [mA]

% - pre-allocating memory
I = zeros(len);
mean_u = zeros(len);

% IC
j = 1;
l = 1;
k = 1+round(dt/(tf/length(tspan)));
I(1) = sigma.*randn*sqrt(2.*dt./tau); % [mA]

mu_width = round(step./dt)+1; % number of time point in 2s

t=(0+dt):dt:tf;
t = t(1:(end-1));

for ii=t
%     disp(ii)
    width = round(ii./dt)+1;
    
    if width <= (mu_width*j - (j-1))
%         sprintf('k=%d',k)
        I(k) = I(k-1) + (mu(l)-I(k-1)).*dt./tau + sigma.*randn*sqrt(2.*dt./tau); % [nA]
        mean_u(k) = mu(l); % current mean
        k = k + round(dt/(tf/length(tspan)));
        if width==(mu_width*j - (j-1))
            j = j+1;
            l = l+1;
        end
    end
end

%% visualization
figure;
plot(tspan, I(1:length(tspan))); hold on;
plot(tspan, mean_u(1:length(tspan)),'r');

end

