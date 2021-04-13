#helper functions


theme_pm <- function () {
  theme_bw(base_size=12) +
    theme(
      panel.grid=element_line(linetype="dashed", color="light grey", size=0.2),
      axis.ticks.length=unit(-0.25, "cm"),
      axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5),
                                             "cm"),color="black", face="bold"),
      axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5),
                                             "cm"),color="black",face="bold"),
      axis.title.x = element_text(face="bold"),
      axis.title.y = element_text(face="bold"),
      title = element_text(face="bold"),
      text=element_text(family = "sans") #sans means arial font for windows users
    )
}
#From http://www.sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization#barplot-with-error-bars
#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summarized
# groupnames : vector of column names to be used as
# grouping variables 
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}








####Run the RunIdyno.py from the command prompt
#inputfile is the path of the xmlfile to be run. Default number of multiples=5
#Changing the order of the runpytho inputs so as to accommodate parallel runs.
runpytho1= function(inputfile,xmlfolder,multiples){
  #the scripts_to_start_idynomics folder has to be the working directory
  python_folder=".../iDynoMiCS/scripts_to_start_idynomics" 
  setwd(python_folder)
  l1=paste("python RunIdyno.py ",xmlfolder,inputfile," --multiples=",multiples,sep="")
  system(l1)
}


####Run the RunIdyno.py from the command prompt
#inputfile is the path of the xmlfile to be run. Default number of multiples=5
#Changing the order of the runpytho inputs so as to accommodate parallel runs.
#for the surfactants folder use iDynoMiCS-master
runpytho_surfactants= function(inputfile,xmlfolder,multiples){
  #the scripts_to_start_idynomics folder has to be the working directory
  python_folder=".../iDynoMiCS-master/scripts_to_start_idynomics" 
  setwd(python_folder)
  l1=paste("python RunIdyno.py ",xmlfolder,inputfile," --multiples=",multiples,sep="")
  system(l1)
}

#To get the total biovolume in the biofilm. 
agent_getTotalBiovolume=function(resultFileFolder,resultFileType,numTimepoints,outputPeriod,columnName){
  unzip('agent_Sum.zip',exdir = paste0(getwd(),'/agent_Sum'))
  species_names=colnames(getSpeciesSpecificAbundance(resultFileFolder,numTimepoints,outputPeriod))
  Totalbiovolume=0
  for (k in 1:length(species_names)) {
    volume=agent_getVolumeOverTime(resultFileFolder,resultFileType,numTimepoints,outputPeriod
                                   ,species_names[k],columnName)
    #If the EPS does not exist in the beginning, the result of agent_getVolumeOverTime will be a single number,
    #that is why we check to see if the number of rows of volume is greater than 1 o.w. there will be a 
    #subscript out of bounds error.
    if(nrow(volume)>1){
      Totalbiovolume=Totalbiovolume+volume[2,]
    }else{
      Totalbiovolume=Totalbiovolume+volume
    }
    
  }
  return(Totalbiovolume)
}

agent_getTotalBiovolume_chemotaxis=function(resultFileFolder,resultFileType,numTimepoints,outputPeriod,columnName){
  unzip('agent_Sum.zip',exdir = paste0(getwd(),'/agent_Sum'))
  species_names=colnames(getSpeciesSpecificAbundance(resultFileFolder,numTimepoints,outputPeriod))
  Totalbiovolume=0
  for (k in 1:length(species_names)) {
    #in the chemotaxis model there are planktonic cells which we wish to avoid. 
    #Species names should start with My, like Mystrep
    if(grepl("My",species_names[k])){
      volume=agent_getVolumeOverTime(resultFileFolder,resultFileType,numTimepoints,outputPeriod
                                     ,species_names[k],columnName)
      #If the EPS does not exist in the beginning, the result of agent_getVolumeOverTime will be a single number,
      #that is why we check to see if the number of rows of volume is greater than 1 o.w. there will be a 
      #subscript out of bounds error.
      if(nrow(volume)>1){
        Totalbiovolume=Totalbiovolume+volume[2,]
      }else{
        Totalbiovolume=Totalbiovolume+volume
      }
      
    }
  }
  
  return(Totalbiovolume)
}



#To get the biovolume of a single species over time

agent_getVolumeOverTime <- function(resultFileFolder,resultFileType,numTimepoints,outputPeriod,speciesReqd,columnName)
{
  allMeasureResults<-NULL
  
  for(tp in seq(from=0,to=numTimepoints,by=outputPeriod))
  {
    # Read in the result
    xmlResultFile<-readSimResultFile(resultFileFolder,resultFileType,tp)
    
    allSpecies<-agent_returnSpeciesResultData(xmlResultFile)
    
    # Now get the measure of the species
    measure<-agent_returnSpeciesVolumeTotal(allSpecies,speciesReqd,columnName)
    
    # Append to the list of growth rates
    allMeasureResults<-rbind(allMeasureResults,measure[[1]])
    
  }
  return(allMeasureResults)
}

#to get the total biovolume for a particular species from a particular output file 
agent_returnSpeciesVolumeTotal =  function(allSpecies,speciesReqd, columnName)
{
  # Get the results for this species name, if there
  speciesResults = data.frame(allSpecies[[speciesReqd]],check.rows = FALSE,row.names=NULL)
  volume=0
  if(nrow(speciesResults)>0)
  {
    for (i in 1:nrow(speciesResults)) {
      volume=volume+(4/3*pi*(speciesResults[i,columnName]^3))
    }
    return(volume)
  }
  else
  {
    return(NULL)
  }
}


#To get the surface area of a single species over time
agent_getsurfaceAreaOverTime <- function(resultFileFolder,resultFileType,numTimepoints,outputPeriod,speciesReqd,columnName)
{
  allMeasureResults<-NULL
  
  for(tp in seq(from=0,to=numTimepoints,by=outputPeriod))
  {
    # Read in the result
    xmlResultFile<-readSimResultFile(resultFileFolder,resultFileType,tp)
    
    allSpecies<-agent_returnSpeciesResultData(xmlResultFile)
    
    # Now get the measure of the species
    measure<-agent_returnSpeciessurfaceAreaTotal(allSpecies,speciesReqd,columnName)
    
    # Append to the list of growth rates
    allMeasureResults<-rbind(allMeasureResults,measure[[1]])
    
  }
  return(allMeasureResults)
}

#to get the total surface area for a particular species from a particular output file 
agent_returnSpeciessurfaceAreaTotal =  function(allSpecies,speciesReqd, columnName)
{
  # Get the results for this species name, if there
  speciesResults = data.frame(allSpecies[[speciesReqd]],check.rows = FALSE,row.names=NULL)
  surfaceArea=0
  if(nrow(speciesResults)>0)
  {
    for (i in 1:nrow(speciesResults)) {
      surfaceArea=surfaceArea+(4*pi*(speciesResults[[columnName]]^2))
    }
    return(surfaceArea)
  }
  else
  {
    return(NULL)
  }
}


#To get the ratio of total surface area to total biovolume in the biofilm. 
agent_getTotalSAtoBV=function(resultFileFolder,resultFileType,numTimepoints,outputPeriod,columnName){
  unzip('agent_Sum.zip',exdir = paste0(getwd(),'/agent_Sum'))
  species_names=colnames(getSpeciesSpecificAbundance(resultFileFolder,numTimepoints,outputPeriod))
  Totalbiovolume=0
  TotalsurfaceArea=0
  for (k in 1:length(species_names)) {
    volume=agent_getVolumeOverTime(resultFileFolder,resultFileType,numTimepoints,outputPeriod
                                   ,species_names[k],columnName)
    surfaceArea=agent_getsurfaceAreaOverTime(resultFileFolder,resultFileType,numTimepoints,outputPeriod
                                             ,species_names[k],columnName)
    #If the EPS does not exist in the beginning, the result of agent_getVolumeOverTime will be a single number,
    #that is why we check to see if the number of rows of volume is greater than 1 o.w. there will be a 
    #subscript out of bounds error.
    if(nrow(volume)>1){
      Totalbiovolume=Totalbiovolume+volume[2,]
      TotalsurfaceArea=TotalsurfaceArea+surfaceArea[2,]
    }else{
      Totalbiovolume=Totalbiovolume+volume
      TotalsurfaceArea=TotalsurfaceArea+surfaceArea
    }
    SAtoBV=TotalsurfaceArea/Totalbiovolume
  }
  return(SAtoBV)
}




speciesnumber_noEPS = function(resultFileFolder,numTimepoints,outputPeriod){
 
  
  speciesAbundance<-getSpeciesSpecificAbundance(resultFileFolder,numTimepoints,outputPeriod)
  
  #to remove the EPS species from the dataframe
  if(any(grepl('EPS',colnames(speciesAbundance)))=='TRUE'){
    speciesAbud=as.data.frame(speciesAbundance[,-which(grepl('EPS',colnames(speciesAbundance)))])
    colnames(speciesAbud)=colnames(speciesAbundance)[-which(grepl('EPS',colnames(speciesAbundance)))]
  }
  else{
    speciesAbud=speciesAbundance
    colnames(speciesAbud)=colnames(speciesAbundance)
  }
  
  return(speciesAbud)
}



