This folder contains the model codes with Ih blocked (noIh) and with Ih present (Ih).

To solve each case, please run:
- [mainIh.m](mainIh.m) for Ih present, and
- [mainNoIh.m](mainNoIh.m) for Ih blocked.

In these scripts, you will find the following stimulations paradigms:
- *Const_current:* when selected -> Somatic input current is a current step pulse, whereas the dendritic input can be either a current step pulse or a double exponential. You can specify the amplitude of both stimulations as well as the duration of the square pulse by defining the stimulation onset and stimulation end. To apply an EPSP-like current or double exponential current to the dendritic compartment, set the stimulation end as 0ms. Dendritic stimulation will start at the time set by the stimulation onset. Please read comments on the codes for more details. These stimulation paradigms were used to study the BAC-firing property shown by these neurons.
- *Pulse_train:*  when selected -> the somatic/dendritic input is a train of square pulses. You can specify their duration, amplitude, and frequency. You can also set a time window for the stimulation -> stimulation onset and end. This type of stimulation was used to study the influence of the somatic APs frequency on the dendritic Ca2+ spike generation.   
- *Noisy staircase current:* when none of the above is selected, a noisy staircase current is applied to the somatic/dendritic compartments. These in vivo-like input currents are generated using the Ornstein-Uhlenbeck method. You can specify its mean, standard deviation, step duration, and step length. You can also choose the compartment in which you want to apply the current. This type of stimulation was used to generate the f-I curve for either somatic or dendritic stimulations.

All parameters in the codes are the ones used to performe the simulations in the present study. Simulation's results are stored in a folder named FiguresData inside the folder called GenerateFigures.
