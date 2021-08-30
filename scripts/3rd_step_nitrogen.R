#This is code for applying the GAM model built in the 1st step on the transects outside
#the Netherlands (non-dutch).The aim is predicting the optimal nitrogen value of the 
#species in each non-dutch transect.
# We then produce histogram to observe the distribution of the predicted values.
#Part of the script includes the old version of the method as developed in the MSc thesis
#Code developed by Alessandro Mari in the period April-October 2020.
#R version 4.0.0


library(here)
library(readxl)
library(raster)
library(ggplot2)

#sourcing the script for building the GAM model
source(here("scripts","1st_step_nitrogen.R"))

#### optimal nitrogen value prediction for non-dutch butterfly species ####

#dataframe where the CNI and average PCs value will be saved for each non-dutch transect 
#year combination
non.dutch.avg <- data.frame(Species_name_2019=character(),
                            bms_id=character(),
                            transect_id=character(),
                            year=numeric(),
                            count=integer(),
                            cni=numeric(),
                            pca1=numeric(),
                            pca2=numeric(),
                            pca3=numeric(),
                            pca4=numeric(),
                            dist_km=numeric(), #distance from centroid 
                            buffer_km = numeric()) 

#filtering all.species.ok with only the species for which we have traits values
all.species.t2 <-  merge(all.species.ok, traits, by.x = "Species_name_2019",
                        by.y = "Species")

#function for calculating the distance between a transect and the NL centroid
nl.dist <- function(lon, lat){
  pointDistance(c(st_coordinates(nl.centr)[1], st_coordinates(nl.centr)[2]), 
                c(lon, lat), lonlat = FALSE)
}


buffer.start <- 160000 #starting buffer size in meters
buffer.end <- 1600000 #ending buffer size in meters
increase <- 60000 #increase of buffer size at each loop repetition
i_max <- (buffer.end - buffer.start)/increase #maximum index for the main loop
index <- 0 #initializing the index for non.dutch.avg rows
#list that will contain  the indexes of transects inside a buffer
inside.list <- vector(mode = "list", length = i_max + 1)

#starting the loop for increasing buffer size at each iteration

for(i in 0:i_max) {
  #increasing the buffer size increases of "increase" after every iteration
  buffer.size <- buffer.start + i * increase 
  buffer <- st_buffer(nl.centr, dist=buffer.size) #buffer object
  lon <- site.locations.ok$transect_lon #transects longitude
  lat <- site.locations.ok$transect_lat  #transects latitude
  buffer.coords <- st_coordinates(buffer) #buffer coordinates
  buffer.lon <- buffer.coords[,1] #buffer longitude points
  buffer.lat <- buffer.coords[,2] #buffer latitude points
  
  #retrieving indexes of the transects between two consecutive buffers (ring)
  
  #index of transects inside the buffer area
  inside.list[[i+1]] <- which(point.in.polygon(lon,lat,buffer.lon,buffer.lat) != 0) 
  if(i == 0){inside <- inside.list[[i+1]]} else {
    inside <- inside.list[[i+1]][!(inside.list[[i+1]] %in% inside.list[[i]])]}
  if(length(inside) > 0){
    #filtering site.locations.ok for transects inside buffer area
    site.locations.in <- site.locations.ok[inside, ] #transects inside buffer area
    #filtering all.species.t2 for the transects inside the buffer area
    all.species.in <- all.species.t2[transect_id %in% site.locations.in$transect_id] #count database filtered by
    #only transects inside the buffer
    
    #starting the loop on bms_id unique values
    
    BMS_unique <- unique(all.species.in$bms_id) #unique values for BMS_id
    for(z in 1:length(BMS_unique)) {
      #filtering all.species.in for a single BMS
      data1 <- all.species.in[bms_id==BMS_unique[z],]
      
      #starting the loop on transect_id unique values
      
      transect_unique <- unique(data1$transect_id)
      for (j in 1:length(transect_unique)) {
        #filtering data1 for a single BMs-transect combination
        data2 <- data1[transect_id==transect_unique[j],]
        t.lon <- site.locations.in[transect_id == transect_unique[j],]$transect_lon
        t.lat <- site.locations.in[transect_id == transect_unique[j],]$transect_lat
        dist <- nl.dist(t.lon, t.lat) #distance of the transect from centroid
        
        #starting the loop on year unique values
        
        year_unique <- unique(data2$year)
        for (k in 1:length(year_unique)) {
          #filtering data1 for a single BMs-transect-year combination
          data3 <- data2[year == year_unique[k],]
         
          #starting the loop on Species_name_2019 looking for the species not included 
          #in nitro.2 and therefore whose optimal nitrogen value is not available
          for(l in 1:length(data3$Species_name_2019)) {
            #allowing entering the loop when the species is not found in nitro.2
            if(all(data3$Species_name_2019[l] != nitro.2$srt_Scientific_name) == "TRUE") {
              
              #old version of the code to filter data3 for species with SSI>11
              #specialist <- data3[data3$SSI > 11, ]
              
              #filtering nitro.train for the species in the transect-year combination
              nitro.val <- nitro.2[srt_Scientific_name %in% data3$Species_name_2019]
              #filtering data3 for species in nitro.val 
              spec.val <- data3[Species_name_2019 %in% nitro.val$srt_Scientific_name]
             
              #calculating average values only if we have more than 1 species in the transect
              if(dim(nitro.val)[1] > 1){
                index <- index + 1 #updating the index for non.dutch.avg rows
                
                #calculating the CNI and the average PCs value in the transect-year combination
                cni <- weighted.mean(nitro.val$Nopt, spec.val$count, na.rm = T)
                pca1 <- weighted.mean(data3$PCA1, data3$count, na.rm = T)
                pca2 <- weighted.mean(data3$PCA2, data3$count, na.rm = T)
                pca3 <- weighted.mean(data3$PCA3, data3$count, na.rm = T)
                pca4 <- weighted.mean(data3$PCA4, data3$count, na.rm = T)
                
                #saving the values in non.dutch.avg
                non.dutch.avg[index,] <- list(data3$Species_name_2019[l], 
                                                BMS_unique[z], 
                                                transect_unique[j], 
                                                year_unique[k], 
                                                data3$count[l],
                                                cni, 
                                                pca1, 
                                                pca2, 
                                                pca3, 
                                                pca4, 
                                                round(dist/1000,0), 
                                                round(buffer.size/1000,0))
              }
              
            }
            
          }
          
        }
      }
    }
  }
}
#converting non.dutch.avg to datatable format
non.dutch.avg <- as.data.table(non.dutch.avg)
#working with a copy of non.dutch.avg
non.dutch.avg2 <- non.dutch.avg
non.dutch.avg2 <- non.dutch.avg2[order(Species_name_2019),]

#averaging CNI and PCs over all the transect/year combination to have a single value
#for each non-dutch species
non.dutch.avg2.sp <-
  non.dutch.avg2 %>% group_by(Species_name_2019) %>% summarize(
    cni_avg = weighted.mean(cni, count, na.rm = T),
    pca1 = weighted.mean(pca1, count, na.rm = T),
    pca2 = weighted.mean(pca2, count, na.rm = T),
    pca3 = weighted.mean(pca3, count, na.rm = T),
    pca4 = weighted.mean(pca4, count, na.rm = T)
  )

#predicting the optimal nitrogen value and standard error using the GAM model "modgam_ok"
N.pred2 <- modgam_ok %>% predict(non.dutch.avg2.sp[,c(2,3,5)], se.fit=T)


#dataframe with predictions and the "observed" optimal nitrogen value for Dutch species
non.dutch.pred <- data.frame(Species_name_2019 = non.dutch.avg2.sp$Species_name_2019,
                             Nopt.gam = round(N.pred2$fit,1),
                             se = round(N.pred2$se.fit,1)) %>% as.data.table()

#### histogram for the optimal nitrogen values of the non-dutch species ####


hist(non.dutch.pred$Nopt.gam, xlab="optimal nitrogen value", ylab="Frequency", col = "grey",
     main = '', cex.lab = 1.3)
abline(v = mean(non.dutch.pred$Nopt.gam), col = "red")
