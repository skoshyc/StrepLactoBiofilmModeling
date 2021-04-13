# StrepLactoBiofilmModeling

All the .xml files are protocol files for iDynoMiCS v 1.2. 
The files with surfactant in the name should be run with the iDynoMiCS version found in https://github.com/alexwweston/iDynoMiCS 

The R code files are as follows:

idynomicshelper.R - Contains functions to run the protocol file using command prompt and python, estimate biovolume, and surface area of the agents from the iDynoMiCS v 1.2 xml result files. Also contains code for a theme for ggplot2 and a function to estimate mean and standard deviation of repeats.

idynomics_solute_plots.R - Code to plot the Toxin/surfactant concentration vs time in a dual biofilm using functions from the iDynoR package.

model_simulations.R - Code to run multiple simulations of the iDynoMiCS models, find mean, standard deviation, and coefficient of variation of the repeats and do a barplot with error bars. 

model_simulations_surfactants.R - Code to run multiple simulations of the iDynoMiCS surfactant models, find mean, standard deviation, and coefficient of variation of the repeats and do a barplot with error bars. 
