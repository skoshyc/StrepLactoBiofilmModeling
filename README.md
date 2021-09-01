# StrepLactoBiofilmModeling

All the .xml files are protocol files for iDynoMiCS v 1.2. 
The files with surfactant in the name should be run with the iDynoMiCS version found in https://github.com/alexwweston/iDynoMiCS 

The R code files are as follows:

idynomicshelper.R - Contains functions to run the protocol file using command prompt and python, estimate biovolume, and surface area of the agents from the iDynoMiCS v 1.2 xml result files. Also contains code for a theme for ggplot2 and a function to estimate mean and standard deviation of repeats.

idynomics_solute_plots.R - Code to plot the Toxin/surfactant concentration vs time in a dual biofilm using functions from the iDynoR package.

model_simulations.R - Code to run multiple simulations of the iDynoMiCS models, find mean, standard deviation, and coefficient of variation of the repeats and do a barplot with error bars. 

model_simulations_surfactants.R - Code to run multiple simulations of the iDynoMiCS surfactant models, find mean, standard deviation, and coefficient of variation of the repeats and do a barplot with error bars. 

The .xml files are as follows:
lacto_single.xml: Protocol file to generate a biofilm with only Lactobacillus paracasei

lacto_single_toxin_production.xml: Protocol file to generate a Lactobacillus paracasei biofilm with production of inhibitory toxin

lacto_single_surfactant_production.xml: Protocol file to generate a Lactobacillus paracasei biofilm with production of surfactant

lacto_single_toxinAndsurfactant_production.xml: Protocol file to generate a Lactobacillus paracasei biofilm with production of inhibitory toxin and a surfactant

strep_single.xml: Protocol file to generate a biofilm with only Streptococcus oralis

strep_lacto_nullmodel.xml: Protocol file to generate S. oralis-L. paracasei biofilm where they compete only for space.

strep_lacto_competition.xml: Protocol file to generate S. oralis-L. paracasei biofilm where they compete for space and nutrients. 

strep_lacto_inhibition.xml: Protocol file to generate S. oralis-L. paracasei biofilm where they compete for space and nutrients and when L. paracasei secretes an inhibitor to the growth of S. oralis. 

strep_lacto_surfactant.xml: Protocol file to generate S. oralis-L. paracasei biofilm where they compete for space and nutrients and when L. paracasei secretes a surfactant. The surfactant causes the biofilm cells of both species to become planktonic.  

strep_lacto_inhibitionAndsurfactant.xml: Protocol file to generate S. oralis-L. paracasei biofilm where they compete for space and nutrients and when L. paracasei secretes a surfactant and an inhibitory substance. The surfactant causes the biofilm cells of both species to become planktonic. The inhibitory substance retards the growth of S. oralis.

