%% test gating variables and time constants

clear
close all
clc

V=-150:0.1:120; % Voltage

%% Calcium L-type current
[CaLmoo,CaLtaum] = CaLKineticsW(V);

%% Na Hodgkin-Huxley type
NaHHmoo = zeros(1, length(V));
NaHHtaum = zeros(1, length(V));
NaHHhoo = zeros(1, length(V));
NaHHtauh = zeros(1, length(V));

for i=1:length(V)
    [NaHHmoo(i),NaHHtaum(i), NaHHhoo(i),NaHHtauh(i)] = NaKinetics(V(i));
end

%% Persistent Na current
Napmoo = zeros(1, length(V));
Naptaum = zeros(1, length(V));
Naphoo = zeros(1, length(V));
Naptauh = zeros(1, length(V));

for i=1:length(V)
    [Napmoo(i),Naptaum(i), Naphoo(i),Naptauh(i)] = NapKinetics(V(i));
end

%% K Hodgkin-Huxley type
KHHmoo = zeros(1, length(V));
KHHtaum = zeros(1, length(V));

for i=1:length(V)
    [KHHmoo(i),KHHtaum(i)] = KKinetics(V(i));
end

%% Slow inactivating K+ current
Kslowmoo = zeros(1, length(V));
Kslowtaum = zeros(1, length(V));
Kslowhoo = zeros(1, length(V));
Kslowtauh = zeros(1, length(V));

for i=1:length(V)
    [Kslowmoo(i),Kslowtaum(i), Kslowhoo(i),Kslowtauh(i)] = KslowKinetics(V(i));
end

%% Ih, Hyperporalized Non-specific cation current
Ihmoo = zeros(1, length(V));
Ihtaum = zeros(1, length(V));

for i=1:length(V)
    [Ihmoo(i),Ihtaum(i)] = IhKinetics(V(i));
end

%% Muscarinic K+ current, Im
[Immoo,Imtaum]=ImKinetics(V);

%% Save data

% save('Kinetics_Correct.mat');

%% Visualization

% Gating variables

figure;

subplot(3, 3, 1)
plot(V,CaLmoo,'r');
legend('m_{\infty}');
title('Ca-Ltype');
xlabel('Voltage (mV)');

subplot(3, 3, 2)
plot(V,NaHHmoo,'r'); hold on; plot(V,NaHHhoo,'b');
legend('m_{\infty}', 'h_{\infty}');
title('Na-HH');
xlabel('Voltage (mV)');

subplot(3, 3, 3)
plot(V,Napmoo,'r'); hold on; plot(V,Naphoo,'b');
legend('m_{\infty}', 'h_{\infty}');
title('Nap');
xlabel('Voltage (mV)');

subplot(3, 3, 4)
plot(V,KHHmoo,'r');
legend('m_{\infty}');
title('K-HH');
xlabel('Voltage (mV)');

subplot(3, 3, 5)
plot(V,Kslowmoo,'r'); hold on; plot(V,Kslowhoo,'b');
legend('m_{\infty}', 'h_{\infty}');
title('Kslow');
xlabel('Voltage (mV)');

subplot(3, 3, 6)
plot(V,Ihmoo,'r');
legend('m_{\infty}');
title('Ih');
xlabel('Voltage (mV)');

subplot(3, 3, 7)
plot(V,Immoo,'r');
legend('m_{\infty}');
title('Im');
xlabel('Voltage (mV)');

%% Time constants

figure;

subplot(3, 3, 1)
plot(V,CaLtaum,'r');
legend('\tau_{m}');
title('Ca-Ltype');
xlabel('Voltage (mV)');
ylabel('seconds');

subplot(3, 3, 2)
plot(V,NaHHtaum,'r'); hold on; plot(V,NaHHtauh,'b');
legend('\tau_{m}', '\tau_{h}');
title('Na-HH');
xlabel('Voltage (mV)');
ylabel('seconds');

subplot(3, 3, 3)
plot(V,Naptaum,'r'); hold on; plot(V,Naptauh,'b');
legend('\tau_{m}', '\tau_{h}');
title('Nap');
xlabel('Voltage (mV)');
ylabel('seconds');

subplot(3, 3, 4)
plot(V,KHHtaum,'r');
legend('\tau_{m}');
title('K-HH');
xlabel('Voltage (mV)');
ylabel('seconds');

subplot(3, 3, 5)
plot(V,Kslowtaum,'r'); hold on; plot(V,Kslowtauh,'b');
legend('\tau_{m}', '\tau_{h}');
title('Kslow');
xlabel('Voltage (mV)');
ylabel('seconds');

subplot(3, 3, 6)
plot(V,Ihtaum,'r');
legend('\tau_{m}');
title('Ih');
xlabel('Voltage (mV)');
ylabel('seconds');

subplot(3, 3, 7)
plot(V,Imtaum.*ones(1,length(V)),'r');
legend('\tau_{m}');
title('Im');
xlabel('Voltage (mV)');
ylabel('seconds');


