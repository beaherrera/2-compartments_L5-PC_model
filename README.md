# 2-compartments_L5-PC_model
Two-compartments Layer 5 Pyramidal Cell Model
- Soma/basal-dendrites compartment includes typical Na+ and K+ conductances. 
- Trunk/apical-dendrites compartments includes: persistent Na+ (Nap), hyperpolarization-activated cation (Ih), slow inactivation K+ (Ks), muscarinic K+ (IM) and Ca2+ L-type currents (CaL). 

## Plataform
MatLab (release 2018b)

A memory of at least 16 GB is needed to run the trials to create the f-I curves. 

## Toolbox
For the numerical simulations and data processing the following MatLab toolboxes were used:
- [CSDplotter-master](https://github.com/espenhgn/CSDplotter)
- [RasterPlot](https://www.mathworks.com/matlabcentral/fileexchange/45671-flexible-and-fast-spike-raster-plotting)
- [SDETools-master](https://github.com/horchler/SDETools)
- [EEGLAB](https://sccn.ucsd.edu/eeglab/index.php), specifically the function [eegfilt.m](https://sccn.ucsd.edu/~arno/eeglab/auto/eegfilt.html).
You can download them from the links or find them in the folder [Toolbox]().

## Description
- [GenerateFigures](GenerateFigures): this folder contains the programs used to create the figures.

## Contributing
Comments and pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
Please make sure to update tests as appropriate.

## License
[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)
