#Running multiple simulations and getting the biovolume/species number from each folder
#This code is specific for the surfactant model and the  inhibition-surfactant combination model
library(vegan)
library(XML)
library(iDynoR)
library(stringr)
library(dplyr)
library(lattice)
library(xml2)
library(gtools)
library(ggplot2)

#always unzip the folders before applying the functions

source(file='.../idynomicshelper.R')

path1="insert path to working directory"


xmlpath_protocol="insert path to protocol file"

xmlpath_results="insert path to result folder"

species="insert name of species for example strep_lacto_surfactant" 

number_species=2#number of species 

xmlpath=paste0(xmlpath_protocol,species,"/")


protocol_file=paste0(species,".xml")

number_repeats=5
#run the specified protocol file number_repeats times
#runpytho_surfactants runs the surfactant and inhibition-surfactant combination mdoel
runpytho_surfactants(protocol_file,xmlpath,number_repeats) 


speciesabud=matrix(0,number_repeats,number_species*2) #number of columns the number of species
#In the surfactant model, for each bacteria there are two type of cells-planktonic and biofilm
biovolume=matrix(0,number_repeats,number_species)
resultfolder=paste0(xmlpath_results,species,"/")
list1=dir(resultfolder);
for(i in 1:number_repeats){
  xmlpath1=paste0(resultfolder,list1[i])
  setwd(xmlpath1)
  unzip('agent_Sum.zip',exdir = paste0(getwd(),'/agent_Sum'))
  unzip('env_State.zip',exdir = paste0(getwd(),'/env_State'))
  unzip('agent_State.zip',exdir = paste0(getwd(),'/agent_State'))
  list2=mixedsort(dir(paste0(xmlpath1,'/agent_State')))
  n=as.numeric(str_extract(list2[length(list2)], "[0-9]+"))
  colnames(speciesabud)=colnames(speciesnumber_noEPS(xmlpath1,n,n))
  if(number_species==1){
    speciesabud[i,] =c(speciesnumber_noEPS(xmlpath1,n,n)[2,]) # to get the species count at last iteration
  }else{
    speciesabud[i,]=as.matrix(speciesnumber_noEPS(xmlpath1,n,n)[2,])
  }
  
  #For biovolume and biomass calculations, we want to take only the bacteria and not the 
  #plank or joined type bacteria
  names <- colnames(speciesabud)[which(grepl('My',colnames(speciesabud)))]
  colnames(biovolume)=paste0(names,"_biovol")
  colnames(biomass)=paste0(names,"_biomass")
  for(k in 1:number_species){
    if(nrow(agent_getVolumeOverTime(xmlpath1,"agent_State",n,n,
                                    names[k],"radius"))>1)
    {biovolume[i,k]=agent_getVolumeOverTime(xmlpath1,"agent_State",n,n,
                                           names[k],"radius")[2,]
    }
    if(speciesabud[i, which(grepl(names[k],colnames(speciesabud)))]!=0 &
       nrow(agent_getVolumeOverTime(xmlpath1,"agent_State",n,n,
                                    names[k],"radius"))==1){
        biovolume[i,k]=agent_getVolumeOverTime(xmlpath1,"agent_State",n,n,
                                               names[k],"radius")
      }
    if(speciesabud[i, which(grepl(names[k],colnames(speciesabud)))]==0){
        biovolume[i,k]=0
      }
  }
  
}


output_matrix=cbind.data.frame(speciesabud,biovolume,"species"=c(rep(species,number_repeats)))
setwd(path1)
write.csv(output_matrix, file = paste0(species,"_surfactants_results_nutrients_only.csv"))





#find the mean+sd of the repetitions
path_folder="insert path to working folder with the files containing biovolume/species count"
setwd(path_folder)
#Create the combined data set
data_lacto=read.csv(file = "insert name of lacto results .csv file",header = TRUE,row.names = 1)
data_strep=read.csv(file = "insert name of strep results .csv file",header = TRUE,row.names = 1)
data_strep_lacto=read.csv(file = "insert name of strep_lacto dual model results .csv file",header = TRUE,row.names = 1)

#Use the relevant type of interaction 
type_interaction="surfactants"
type_interaction="inhibitionAndsurfactants"

filename=paste0("combined_",type_interaction,"_biovolume_")
filename_species_count=paste0("combined_",type_interaction,"_species_count_")

number_repeats=5
species1=c(rep("Lacto",number_repeats))
species2=c(rep("Strep",number_repeats))

combined_dataframe=cbind.data.frame("biovol"=data_lacto$Mylacto_biovol,"species"=species1,
                                    "biofilm"=c(rep("Lacto",number_repeats)))

combined_dataframe=rbind.data.frame(combined_dataframe,
                                    cbind.data.frame("biovol"=data_strep$Mystrep_biovol,"species"=species2,
                                                     "biofilm"=c(rep("Strep",number_repeats))))
#change mixed species biofilm to the appropriate dataset
combined_dataframe=rbind.data.frame(combined_dataframe,
                                    cbind.data.frame("biovol"=data_strep_lacto$Mylacto_biovol,"species"=species1,
                                                     "biofilm"=c(rep("Lactodual",number_repeats))))
combined_dataframe=rbind.data.frame(combined_dataframe,
                                    cbind.data.frame("biovol"=data_strep_lacto$Mystrep_biovol,"species"=species2,
                                                     "biofilm"=c(rep("Strepdual",number_repeats))))
combined_dataframe=rbind.data.frame(combined_dataframe,
                                    cbind.data.frame("biovol"=(data_strep_lacto$Mystrep_biovol)+
                                                       data_strep_lacto$Mylacto_biovol,
                                                     "species"=c(rep("total_biovol",number_repeats)),
                                                     "biofilm"=c(rep("dual",number_repeats))))
write.csv(combined_dataframe,file=paste0(filename,number_repeats,"runs.csv"))


#number of species
combined_dataframe_species=cbind.data.frame("species_count"=data_lacto$Mylacto,"species"=species1,
                                            "biofilm"=c(rep("Lacto",number_repeats)))

combined_dataframe_species=rbind.data.frame(combined_dataframe_species,
                                            cbind.data.frame("species_count"=data_strep$Mystrep,"species"=species2,
                                                             "biofilm"=c(rep("Strep",number_repeats))))
#change mixed species biofilm to the appropriate dataset
combined_dataframe_species=rbind.data.frame(combined_dataframe_species,
                                            cbind.data.frame("species_count"=data_strep_lacto$Mylacto,"species"=species1,
                                                             "biofilm"=c(rep("Lactodual",number_repeats))))
combined_dataframe_species=rbind.data.frame(combined_dataframe_species,
                                            cbind.data.frame("species_count"=data_strep_lacto$Mystrep,"species"=species2,
                                                             "biofilm"=c(rep("Strepdual",number_repeats))))
combined_dataframe_species=rbind.data.frame(combined_dataframe_species,
                                            cbind.data.frame("species_count"=(data_strep_lacto$Mystrep)+
                                                               data_strep_lacto$Mylacto,
                                                             "species"=c(rep("total_count",number_repeats)),
                                                             "biofilm"=c(rep("dual",number_repeats))))
write.csv(combined_dataframe_species,file=paste0(filename_species_count,number_repeats,"runs.csv"))

#Read in the data and generate the error bar plots
setwd(path_folder)
combined_data=read.csv(file=paste0(filename,number_repeats,"runs.csv"),header=TRUE,row.names = 1)
combined_data_mean= data_summary(combined_data, varname="biovol", 
                                 groupnames=c("species", "biofilm"))
#combined_data_CV includes the coefficient of variation as well
combined_data_CV=cbind.data.frame(combined_data_mean,
                                  "coeff_variation"=combined_data_mean$sd/combined_data_mean$biovol)
write.csv(combined_data_CV,file=paste0(filename,"CV_",number_repeats,"runs.csv"))
# Convert biofilm to a factor variable for the plot
combined_data_mean$biofilm=as.factor(combined_data_mean$biofilm)
# Default bar plot
pdf(paste0(filename,number_repeats,"runs.pdf"))
p = ggplot(combined_data_mean, aes(x=biofilm, y=biovol, fill=species)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=biovol-sd, ymax=biovol+sd), width=.2,
                position=position_dodge(.9))+theme_pm()+
  geom_text(aes(label=species),vjust=.5,hjust=2,angle=90, color="black",
            position = position_dodge(0.9), size=2.5)
p+labs(title="Biovolume for single and double biofilms", y="Biovolume (micron^3)", x = "Type of biofilm")
dev.off()

#species count
combined_data_species=read.csv(file=paste0(filename_species_count,number_repeats,"runs.csv"),header=TRUE,row.names = 1)
combined_data_species_mean= data_summary(combined_data_species, varname="species_count", 
                                         groupnames=c("species", "biofilm"))
#combined_data_CV includes the coefficient of variation as well
combined_data_species_CV=cbind.data.frame(combined_data_species_mean,
                                          "coeff_variation"=combined_data_species_mean$sd/combined_data_species_mean$species_count)
write.csv(combined_data_species_CV,file=paste0(filename_species_count,"CV_",number_repeats,"runs.csv"))
# Convert biofilm to a factor variable for the plot
combined_data_species_mean$biofilm=as.factor(combined_data_species_mean$biofilm)
# Default bar plot
pdf(paste0(filename_species_count,number_repeats,"runs.pdf"))
p = ggplot(combined_data_species_mean, aes(x=biofilm, y=species_count, fill=species)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=(species_count)-sd, ymax=(species_count)+sd), width=.2,
                position=position_dodge(.9))+theme_pm()+
  geom_text(aes(label=species),vjust=.5,hjust=2,angle=90, color="black",
            position = position_dodge(0.9), size=2.5)
p+labs(title="Species count for single and double biofilms", y="Species count", x = "Type of biofilm")
dev.off()

