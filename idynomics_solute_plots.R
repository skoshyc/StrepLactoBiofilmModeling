#Plots the solute concentrations at specified locations across the simulations
library(vegan)
library(XML)
library(iDynoR)
library(stringr)
library(dplyr)
library(lattice)
library(xml2)
library(gtools)
library(ggplot2)
library(reshape2)
#always unzip the folders before applying the functions

source(file='.../idynomicshelper.R')

path1="insert path to working directory"
xmlpath_results="insert path to result folder"


species="insert name of species for example strep_lacto_inhibition" 
number_species= 2#number of species 

resultfolder=paste0(xmlpath_results,species,"/")
list1=dir(resultfolder) #Get the names of all the different repeats
resultfileFolder=file.path(resultfolder,list1[1]) #We took the results from the first folder, thus
#list1[1], it can be changed to desired number
list2=mixedsort(dir(paste0(resultfileFolder,'/env_State')))
#now will get the solute concentration from the locations upper-left:(1,1),upper-right:(1,ncol),
#lower-left: (nrow,1), lower-right:(nrow,ncol), and
#middle:(floor(nrow/2),floor(ncol/2)),  
output_matrix=matrix(0,length(list2),6,dimnames = list(c(),c("Time_h","Upper_left","Upper_right",
                                                             "Lower_left","Lower_right",
                                                             "Middle")))
for (sim_number in 1:length(list2)) {
  timepoint=as.numeric(str_extract(list2[sim_number], "[0-9]+"))
  #plot the solutes at specified locations
  xmlResultFile<-readSimResultFile(resultfileFolder,"env_State",timepoint)
  time_in_hours <- agent_returnSimTime(xmlResultFile)
  soluteReqd=4 #the surfactant/toxin is the fourth soulte specified in the protocol file when it is 
  #only surfactant or toxin produced
  soluteData<-env_returnSpecifiedSoluteData(xmlResultFile,soluteReqd)
  #from the plotContour function of iDynoR
  nums<-c(1:nrow(soluteData))
  minix <- min(soluteData[,2])
  maxix <- max(soluteData[,2])
  x <- minix:maxix 
  miniy <- min(soluteData[,3])
  maxiy <- max(soluteData[,3])
  y <- miniy:maxiy
  temp <- matrix(nrow=length(y), ncol=length(x))
  for(i in x)
  {
    for(j in y)
    {
      if(minix==0)
      {
        iv <- i+1
      }
      else
      {
        iv <- i
      }
      
      if(miniy==0)
      {
        jv <- j+1
      }
      else
      {
        jv <- j
      }
      
      temp[iv,jv] <- soluteData[,5][which(soluteData[,2]==i)][which(soluteData[,3][which(soluteData[,2]==i)]==j)]
    }
  }
  
  
  #instead of heatplot in plotContour, saving results from pre-defined locations
  #the selected locations in the grid are upper-left:(1,1),upper-right:(1,ncol),
  #lower-left: (nrow,1), lower-right:(nrow,ncol), and
  #middle:(floor(nrow/2),floor(ncol/2))
  
  output_matrix[sim_number,1]=time_in_hours
  output_matrix[sim_number,2:6]=c(temp[1,1],temp[1,ncol(temp)],
                                  temp[nrow(temp),1],temp[nrow(temp),ncol(temp)],
                                   temp[floor(nrow(temp)/2),floor(ncol(temp)/2)])
}


setwd(path1)
write.csv(output_matrix,file="Insert desired file name.csv")

df <- as.data.frame(output_matrix)
df <- melt(df ,  id.vars = 'Time_h', variable.name = 'Grid_locations')


#Plot for Toxin and K_I in the inhibition model
png(file=paste0("Plot_",species,"_Toxin_concentration_1stFolder.png"),width = 590, height = 350)
ggplot(df[which(df[,"value"]!=0),], aes(Time_h,value)) + geom_line(aes(colour = Grid_locations))+
  geom_hline(yintercept=0.0025, linetype="dashed", color = "red")+ theme_pm()+
#geom_hline is for the value of K_I=0.0025 for toxin in the inhibition model
  labs(title="Inhibitor concentration vs time", y="Inhibitor concentration (g/L)",x="Time(h)")+
  annotate("text", min(df$Time_h), 0.0025, vjust=-1,label = "K_I")
dev.off()

#Plot for the minimum threshold in the surfactant model
png(file=paste0("Plot_",species_surfactant,"_surfactant_concentration.png"),width = 590, height = 350)

# plot on same grid, each series colored differently -- 
# good if the series have same scale
ggplot(df[which(df[,"value"]!=0),], aes(Time_h,value)) + geom_line(aes(colour = Grid_locations))+
  geom_hline(yintercept=0.005, linetype="dashed", color = "red")+ 
  theme_pm()+
  #geom_hline is for the surfactant model's minimum threshold is 0.005
  labs(title="Surfactant concentration vs time", y="Surfactant concentration (g/L)",x="Time(h)")+
  annotate("text", min(df$Time_h), 0.005, vjust=-1, hjust= 0,label = "Min. tolerance")
dev.off()

#plot of species count
#can use plotTimeCourseAbund
