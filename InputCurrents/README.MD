This folder contains functions that generate different types of imput currents.

- [Inoise.m](Inoise.m): generates a noisy current that is Gaussian distributed with mean mu, standard deviation sigma, and a correlation length tau
- [Inoise_incMu.m](Inoise_incMu.m): creates an in vivo-like inputs using the Ornstein-Uhlenbeck method.

Noisy currents with mean mu, standard deviation sigma, and correction length tau are generated according to the equation:

 I(t+dt)=I(t)+(mu-I(t))/tau+sigma *Gt *sqrt(2*dt/tau)
 
where: Gt is a random number taken each time step from a Gaussian distribution with mean and standard deviation 1.

- [doubleExpFunct.m](doubleExpFunct.m): generates a double exponential functions with time constants tau1 and tau2.
- [square_pulse.m](square_pulse.m): creates a train of square pulses.
