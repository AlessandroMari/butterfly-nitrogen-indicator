#This is code to cross-validate the GAM model built in the 1st step on Dutch transects and
#species.
#Therefore, we predict the optimal nitrogen value for each Dutch species for which
#we have a value already available, and then comparable, in nitro.2. We then produce
#a plot for comparing the linear regression line of predicted vs observed values to the
#1:1 line, which would represent a perfect prediction.
#Part of the script includes the old version of the method as developed in the MSc thesis
#Code developed by Alessandro Mari in the period April-October 2020.
#R version 4.0.0


library(here)
library(readxl)
library(ggplot2)

#sourcing the script for building the GAM model
source(here("scripts","1st_step_nitrogen.R"))


#### Leave-One-Out cross-validation on Dutch species ####

#dataframe where the CNI and average PCs value will be saved for each dutch transect and year
#combination
dutch.avg <- data.frame(Species_name_2019 = character(),
                          bms_id = character(),
                          transect_id = character(),
                          year=numeric(),
                          count = numeric(),
                          cni = numeric(),
                          pca1 = numeric(),
                          pca2 = numeric(),
                          pca3 = numeric(),
                          pca4 = numeric())


index <- 0 #initializing the index for cni_data rows

#starting the loop for the Leave-One-Out Cross_Validation (LOOCV)
for(k in 1:35){ #35 is the number of Dutch species whose Nopt is available
  print(k) #printing k for tracking the loop state
  nitro.test <- nitro.2$srt_Scientific_name[k] #test species whose Nopt gets calculated
  nitro.train <- nitro.2[-k,] #train dataset 
  
  #starting the loop on transect_id unique values
  
  transect_unique <- unique(all.species.t$transect_id)
  for (j in 1:length(transect_unique)) {
    #filtering all.species.in for a single transect
    data1 <- all.species.t[transect_id==transect_unique[j],]
    
    #starting the loop on year unique values
    
    year_unique <- unique(data1$year)
    for (k in 1:length(year_unique)) {
      #filtering the data of all species counts for a single transect-year combination
      data2 <- data1[year == year_unique[k],]
      
      #starting the loop on Species_name_2019 looking for the species not included in the
      #train dataset and therefore equal to nitro.test
      for(l in 1:length(data2$Species_name_2019)) {
        #allowing entering the loop when the species is equal to nitro.test
        if(all(data2$Species_name_2019[l] == nitro.test) == "TRUE") {
          
          #old version of the code to filter data2 for species with SSI>11
          #specialist <- data2[data2$SSI > 11, ]
          
          #filtering nitro.train for the species in the transect-year combination
          nitro.val <- nitro.train[srt_Scientific_name %in% data2$Species_name_2019]
          #filtering data2 for species in nitro.val 
          spec.val <- data2[Species_name_2019 %in% nitro.val$srt_Scientific_name]
          
          #calculating average values only if we have more than 1 species in the transect
          if(dim(nitro.val)[1] > 1){
            index <- index + 1 #updating the index for cni_data rows
            
            #calculating the CNI and the average PCs value in the transect-year combination
            cni <- weighted.mean(nitro.val$Nopt, spec.val$count, na.rm = T)
            pca1 <- weighted.mean(data2$PCA1, data2$count, na.rm = T)
            pca2 <- weighted.mean(data2$PCA2, data2$count, na.rm = T)
            pca3 <- weighted.mean(data2$PCA3, data2$count, na.rm = T)
            pca4 <- weighted.mean(data2$PCA4, data2$count, na.rm = T)
            
            #saving the values in dutch.pred
            dutch.avg[index, ] <- list(data2$Species_name_2019[l], 
                                         "NLBMS",
                                         transect_unique[j],
                                         year_unique[k],
                                         data2$count[l],
                                         cni,
                                         pca1,
                                         pca2,
                                         pca3,
                                         pca4)
            
          }
        }
      }
    }
  }
}
dutch.avg <- as.data.table(dutch.avg) #converting dutch.avg to datatable format

#working with a copy of dutch.avg
dutch.avg.2 <- dutch.avg
dutch.avg.2 <- dutch.avg.2[order(Species_name_2019),]

#averaging CNI and PCs over all the transect/year combination to have a single value
#for each Dutch species
dutch.avg.sp <-
  dutch.avg.2 %>% group_by(Species_name_2019) %>% 
  summarize(
    cni_avg = weighted.mean(cni, count, na.rm = T),
    pca1 = weighted.mean(pca1, count, na.rm =  T),
    pca2 = weighted.mean(pca2, count, na.rm = T),
    pca3 = weighted.mean(pca3, count, na.rm = T),
    pca4 = weighted.mean(pca4, count, na.rm = T)
  )

#predicting the optimal nitrogen value and standard error using the GAM model "modgam_ok"
N.pred <- modgam_ok %>% predict(dutch.avg.sp[,c(2,3,5)], se.fit=T)
#dataframe with predictions and the "observed" optimal nitrogen value for Dutch species
dutch.pred <- data.frame(Species_name_2019 = dutch.avg.sp$Species_name_2019,
                         Nopt.obs = nitro.2$Nopt,
                         Nopt.gam = round(N.pred$fit,1),
                             se = round(N.pred$se.fit,1)) %>% as.data.table()



#### plotting predicted vs observed optimal nitrogen values  ####

mod <- lm(Nopt.obs ~ Nopt.gam , data = dutch.pred)

#prediction plot for the results of modgam1 application

ggplot(dutch.pred,aes(x=Nopt.gam, y=Nopt.obs)) +
  theme_bw() + 
  geom_point() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=18))+
  # geom_smooth(method="lm", color = "red",se=F, fullrange=T) + 
  geom_errorbar(aes(xmin=Nopt.gam-se, xmax=Nopt.gam+se), width=0.2, color = "darkgrey")+
  xlab("predicted optimal nitrogen value")+
  ylab("optimal nitrogen value") +
  xlim(c(0,8)) + 
  ylim(c(0,8)) +
  # theme(axis.text=element_text(size=14),
  #      axis.title=element_text(size=14)) + 
  geom_abline(slope = mod$coefficients[2], intercept = mod$coefficients[1],
              color = "black") +  #regression line in black
  annotate(geom = "text", 
           label = paste("y = ", round(mod$coefficients[2],2),"x ",
                                        format(round(mod$coefficients[1], 2),
                                               nsmall = 2), "\n",
                                        "p-value < 0.001"),
           x = 1, y = 7, size = 7) +
  geom_abline(slope = 1, intercept = 0, color="red") #1:1 line in red



