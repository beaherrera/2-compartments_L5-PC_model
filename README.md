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
- [InputCurrents](InputCurrents): contains the functions used to generate the input currents.
- [IonicCurrents](IonicCurrents): contains the channels kinetics and a test function to visualize the channel kinetics curves: activation/inactivation and time constants.
- [SingleNeuron](SingleNeuron): contains the model codes with Ih blocked and without blocking Ih for a single neuron. **Add [InputCurrents](InputCurrents) and [IonicCurrents](IonicCurrents) to MATLAB path before running scripts.**
- [L5PCs](L5PCs): contains the model codes with Ih blocked and without blocking Ih for a collaction of L5-PCs located randomly in a cortical column. It also contains the codes for the local field potentials (LFPs) and current source density (CSD) estimations. **Add [InputCurrents](InputCurrents) and [IonicCurrents](IonicCurrents) to MATLAB path before running scripts.**
- [GenerateFigures](GenerateFigures): contains the programs used to create the figures.

## License
Copyright (C) 2020 Beatriz Herrera, Amirsaman Sajad, Geoffrey F. Woodman, Jeffrey D. Schall, Jorge J. Riera.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)

## Citation
Herrera B, Sajad A, Woodman GF, Schall JD, Riera JJ. A minimal biophysical model of neocortical pyramidal cells: implications for frontal cortex microcircuitry and field potential generation. J Neurosci **40**, 8513-8529 (2020); DOI: https://doi.org/10.1523/JNEUROSCI.0221-20.2020.
