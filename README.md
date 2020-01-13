# 2-compartments_L5-PC_model
Two-compartments Layer 5 Pyramidal Cell Model

The basal-dendritic/somatic compartment includes the classic Hodgkin-Huxley sodium (I_Na) and potassium delayed rectifier (I_Kdr) currents (Hodgkin and Huxley, 1952). The apical-dendrite/trunk compartment includes persistent Na^+ current (I_Nap) (Magistretti and Alonso, 1999), Ca^(2+) L-type current (I_CaL) (Lytton and Sejnowski, 1991), hyperpolarization-activated non-specific cation current (I_h) (Kole et al., 2006), muscarinic K+ current (I_M) (Adams et al., 1982), and the slow-inactivating potassium current (I_Ks) (Korngreen and Sakmann, 2000).

## Plataform
MatLab (release 2018b)

A memory of at least 16 GB is needed to run the trials to create the f-I curves. 

## Toolbox
For the numerical simulations and data processing the following MatLab toolboxes were used:
- [CSDplotter-master](https://github.com/espenhgn/CSDplotter)
- [RasterPlot](https://www.mathworks.com/matlabcentral/fileexchange/45671-flexible-and-fast-spike-raster-plotting)
- [SDETools-master](https://github.com/horchler/SDETools)
- [EEGLAB](https://sccn.ucsd.edu/eeglab/index.php), specifically the function [eegfilt.m](https://sccn.ucsd.edu/~arno/eeglab/auto/eegfilt.html).
You can download them from the links or find them in the folder [Toolbox](Toolbox).

## Description
- [GenerateFigures](GenerateFigures): contains the programs used to create the figures.
- [InputCurrents](InputCurrents): contains the functions used to generate the input currents.
- [IonicCurrents](IonicCurrents): contains the channels kinetics and a test function to visualize the channel kinetics curves: activation/inactivation and time constants.
- [SingleNeuron](SingleNeuron): contains the model codes with Ih blocked and without blocking Ih for a single neuron. 

## Contributing
Comments and pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
Please make sure to update tests as appropriate.

## License
[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)
