function y = doubleExpFunct(t, tau1, tau2)
% doubleExpFunct -> generates a double exponential functions with time
% constants tau1 and tau2.

y = (1 - exp(-t./tau1)).*exp(-t/tau2);

end