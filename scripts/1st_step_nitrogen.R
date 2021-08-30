#This is code to build a GAM model for predicting butterfly optimal nitrogen value
#based on:
#- optimal nitrogen values available for Dutch species
#- species biological/ecological traits summarized with a PCA
#Part of the script includes the old version of the method as developed in the MSc thesis
#Code developed by Alessandro Mari in the period April-October 2020.
#R version 4.0.0

library(here)
library(readxl)
library(mgcv)

#sourcing the script which invokes data and prepares it
source(here("scripts","butterfly_data.R"))

#reading the csv file with traits data

#old traits data used for the MSc thesis project from Essens et al. (2017)
#th_traits <- read_excel(here("documents", "12258_Essens_traits.xlsx"), 
#                        sheet = "Suppl Table 1", skip = 8)
#th_traits <- as.data.table(th_traits)
#selecting only data relative to species name, SSI and traits principal components 
#th_traits <- th_traits[,c(2,8,10,11,12)]
#names(th_traits)[3:5] <- c("PC_LH1", "PC_LH2", "PC_LH3")
#th_traits <- th_traits[order(Species),]

#new traits data from Wallisdevries (2014)
traits <- fread(here("documents","Butterflies_Traits_Wallisdevries.CSV"))
#selecting only data relative to species name and traits principal components 
traits <- traits[,c(2,20:23)]             


#filtering all.species for only Dutch transects
all.species.2 <- all.species[bms_id == "NLBMS",]
#filtering site.locations for only Dutch transects
site.nl <- site.locations[bms_id == "NLBMS",]

#deleting all the species with a Nopt of NA from nitrod.indicators 
nitro.2 <- nitrod.indicators[!is.na(nitrod.indicators$Nopt),]
#considering only the Dutch species of nitro.2 common to all.species.2 
nitro.2 <- nitro.2[srt_Scientific_name %in% unique(all.species.2$Species_name_2019), ] 
#considering only the Dutch species in all.species.2 common to nitro.2 (for which 
#we have a Nopt value)
#all.species.2 <- all.species.2[Species_name_2019 %in% nitro.2$srt_Scientific_name,]


#merging all.species.2 with th_traits (old version)
#all.species.t <- merge(all.species.2, th_traits, by.x = "Species_name_2019", 
#by.y = "Species")
#excluding the most generalist species (SSI > 11, MSc thesis version)
#all.species.t <- all.species.t[SSI > 11, ]

#merging all.species.2 with traits
all.species.t <- merge(all.species.2, traits, by.x = "Species_name_2019",
                       by.y = "Species")


#calculating the CNI value and the average value of each PC in each transect 
#in the period 2010-2017 

#dataframe where the CNI and average PCs value will be saved for each transect and year
#combination
cni_data <- data.frame(bms_id = character(),
                       transect_id = character(),
                       year=numeric(),
                       cni = numeric(),
                       pca1 = numeric(),
                       pca2 = numeric(),
                       pca3 = numeric(),
                       pca4 = numeric()) 

buffer.start <- 5000 #starting buffer size in meters
buffer.end <- 170000 #ending buffer size in meters
increase <- 165000 #increase of buffer size at each loop repetition
i_max <- (buffer.end - buffer.start)/increase #maximum index for the main loop
index <- 0 #initializing the index for cni_data rows
#list that will contain  the indexes of transects inside a buffer
inside.list <- vector(mode = "list", length = i_max + 1)

#starting the loop for increasing buffer size at each iteration

for(i in 0:i_max) {
  #increasing the buffer size increases of "increase" after every iteration  
  buffer.size <- buffer.start + i * increase 
  buffer <- st_buffer(nl.centr, dist=buffer.size) #buffer object
  lon <- site.nl$transect_lon #transects longitude
  lat <- site.nl$transect_lat  #transects latitude
  buffer.coords <- st_coordinates(buffer) #buffer coordinates
  buffer.lon <- buffer.coords[,1] #buffer longitude points
  buffer.lat <- buffer.coords[,2] #buffer latitude points
  
  #retrieving indexes of the transects between two consecutive buffers (ring)
  
  #index of transects inside the buffer area
  inside.list[[i+1]] <- which(point.in.polygon(lon,lat,buffer.lon,buffer.lat) != 0) 
  if(i == 0){inside <- inside.list[[i+1]]} else {
    inside <- inside.list[[i+1]][!(inside.list[[i+1]] %in% inside.list[[i]])]}
  if(length(inside) > 0){
    #filtering site.nl for transects inside buffer area
    site.locations.in <- site.nl[inside, ] 
    #filtering all.species.t for the transects inside the buffer area
    all.species.in <- all.species.t[transect_id %in% site.locations.in$transect_id] 
    
    #starting the loop on transect_id unique values
   
      transect_unique <- unique(all.species.in$transect_id) 
      for (j in 1:length(transect_unique)) {
        #filtering all.species.in for a single transect
        data1 <- all.species.in[transect_id==transect_unique[j],]
       
        #starting the loop on year unique values
        
        year_unique <- unique(data1$year)
        for (k in 1:length(year_unique)) {
          #filtering data1 for a single transect-year combination
          data2 <- data1[year == year_unique[k],]
          
          index <- index + 1
          
          #filtering nitro.2 for the species in the transect-year combination
          nitro.3 <- nitro.2[srt_Scientific_name %in% data2$Species_name_2019]
          #filtering data2 for species in nitro.3 
          counts <- data2[data2$Species_name_2019 %in% nitro.3$srt_Scientific_name]
          
          #calculating the CNI and the average PCs value in the transect-year combination
          cni <- weighted.mean(nitro.3$Nopt,counts$count, na.rm=T)
          pca1 <- weighted.mean(data2$PCA1, data2$count, na.rm=T)
          pca2 <- weighted.mean(data2$PCA2, data2$count, na.rm=T)
          pca3 <- weighted.mean(data2$PCA3, data2$count, na.rm=T)
          pca4 <- weighted.mean(data2$PCA4, data2$count, na.rm=T)
          
          #saving the values in cni_data
          cni_data[index,] <- list("NLBMS", 
                                   transect_unique[j],
                                   year_unique[k],
                                   round(cni,1),
                                   pca1,
                                   pca2,
                                   pca3,
                                   pca4)
          
        }
    }
  }
}

#converting cni_data to datatable format
cni_data <- as.data.table(cni_data)
#merging all.species.2 with cni_data by each transect-year combination
merge.df <- merge(all.species.2, cni_data, by=c("transect_id","year"))
#merging merge.df with nitro.2 by species name
merge.df <- merge(merge.df, nitro.2[,c(1,2)], by.x= "Species_name_2019", 
                  by.y = "srt_Scientific_name")

#datatable with average values of CNI and (averaged) PCs across all transect-year 
#combination for each Dutch species
cni.avg <-
  merge.df %>% group_by(Species_name_2019, Nopt) %>%  summarize(
    cni_avg = weighted.mean(cni, count, na.rm = T),
    pca1 = weighted.mean(pca1, count, na.rm = T),
    pca2 = weighted.mean(pca2, count, na.rm = T),
    pca3 = weighted.mean(pca3, count, na.rm = T),
    pca4 = weighted.mean(pca4, count, na.rm = T)
  )


#GAM model with Nopt as response variable and cni_avg, pca1 pca2, pca3, and pca4 as
#predictors (explanatory variables)
modgam <- gam(Nopt ~ cni_avg + s(pca1, k=4) + s(pca2, k=4) +s(pca3, k=4) + s(pca4, k=4),
               data  = cni.avg, method = "REML", select=TRUE)
modgam_ok <- gam(Nopt ~ cni_avg + s(pca1, k=4) + s(pca3, k=4),
                 data  = cni.avg, method = "REML", select=TRUE)

#summary of modgam
summary(modgam_ok)

