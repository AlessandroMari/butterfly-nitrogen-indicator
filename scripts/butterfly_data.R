#This is code to import and prepare the data relative to:
#- butterflies occurrences inside the Netherlands
#- butterflies occurrences outside the Netherlands, in Western Europe and filtered for
#  the Atlantic and Continental biogeographical regions
#Code developed by Alessandro Mari in the period April-October 2020.
#R version 4.0.0

#loading required packages
library(here)
library(data.table)
library(sf)
library(rgdal)
library(dplyr)

#datatable with information on butterfly species occurrences in Western Europe
all.species <- readRDS(here("data","Ntot_per_species_per_site_year.rds"))
#ordering all.species for species name
all.species <- all.species[order(Species_name_2019)]

#datatable with spatial information of transects location
site.locations <- readRDS(here("data","site_locations.rds"))

#datatable with the optimal nitrogen value for Dutch butterflies
nitrod.indicators <- fread(here("data","N_indicator_Europe_decimal.csv"), header = T)
#ordering nitrod.indicators for species name
nitrod.indicators <- nitrod.indicators[order(srt_Scientific_name)]

#updating the scientific name of some species
all.species$Species_name_2019 <- ifelse(all.species$Species_name_2019=="Erebia rondoui", "Erebia hispania", all.species$Species_name_2019)
all.species$Species_name_2019 <- ifelse(all.species$Species_name_2019=="Iolana debilitata", "Iolana iolas", all.species$Species_name_2019)
all.species$Species_name_2019 <- ifelse(all.species$Species_name_2019=="Iphiclides feisthamelii", "Iphiclides podalirius", all.species$Species_name_2019)
all.species$Species_name_2019 <- ifelse(all.species$Species_name_2019=="Leptidea reali", "Leptidea sinapis", all.species$Species_name_2019)
all.species$Species_name_2019 <- ifelse(all.species$Species_name_2019=="Leptidea juvernica", "Leptidea sinapis", all.species$Species_name_2019)
all.species$Species_name_2019 <- ifelse(all.species$Species_name_2019=="Lycaena bleusei", "Lycaena tityrus", all.species$Species_name_2019)
all.species$Species_name_2019 <- ifelse(all.species$Species_name_2019=="Polyommatus celina", "Polyommatus icarus", all.species$Species_name_2019)
all.species$Species_name_2019 <- ifelse(all.species$Species_name_2019=="Pontia edusa", "Pontia daplidice", all.species$Species_name_2019)
all.species$Species_name_2019 <- ifelse(all.species$Species_name_2019=="Pyrgus malvoides", "Pyrgus malvae", all.species$Species_name_2019)

#removing migrant Vanessa cardui
all.species <- all.species[Species_name_2019 != "Vanessa cardui"]
#considering only the species that have more than 0 counts
all.species <- all.species[count != 0]
#considering only years from 2010 to 2017
all.species <- all.species[year %in% 2010:2017]

#summing counts for duplicates after changing the old names to new names
all.species <- all.species %>% 
  group_by(Species_name_2019, bms_id, transect_id, year) %>% 
  summarize(count = sum(count)) %>% as.data.table()


#reading the shapefile in EPSG:3035 with European countries 
eu.sf <- st_read(here("data","NUTS_10M_2016_3035_countries.shp"))

#calculating the centroid for the Netherlands
nl.sf <- eu.sf[eu.sf$CNTR_CODE == "NL",]
nl.centr <- st_point_on_surface(nl.sf)

#excluding Dutch transects from all.species and site.locations datatables
all.species.new <- all.species[bms_id != "NLBMS"]
site.locations.new <- site.locations[bms_id != "NLBMS"]

#reading shapefile with EU biogeographical regions
bioreg <- readOGR(here("data","BiogeoRegions2016.shp"))
#transects filtered for Continental and Atlantic biogeographical regions
my.bioreg <- bioreg[bioreg$code %in% c("Atlantic", "Continental"),] 

#vectors with transects longitude and latitude information
lon <- site.locations.new$transect_lon #transects longitude
lat <- site.locations.new$transect_lat  #transects latitude


#retrieving transects inside the Atlantic biogeographical region

#initializing the index b and the vector index.atl
b = 0
index.atl <- vector("numeric", 3721)

for(j in 1:6902){
  polyg.lon <- my.bioreg@polygons[[1]]@Polygons[[j]]@coords[,1]
  polyg.lat <- my.bioreg@polygons[[1]]@Polygons[[j]]@coords[,2]
  index <- which(point.in.polygon(lon,lat,polyg.lon,polyg.lat) != 0)
  
  if(length(index>0)){
    #print(index)
    a <- b + 1 
    b <- a + length(index)  - 1
    index.atl[a:b] <- index 
    
  } 
}


#retrieving transects inside the Continental biogeographical region

#initializing the index b and the vector index.con
b = 0
index.con <- vector("numeric", 3721)

for(j in 1:2592){
  polyg.lon <- my.bioreg@polygons[[2]]@Polygons[[j]]@coords[,1]
  polyg.lat <- my.bioreg@polygons[[2]]@Polygons[[j]]@coords[,2]
  index <- which(point.in.polygon(lon,lat,polyg.lon,polyg.lat) != 0)
  
  if(length(index>0)){
    #print(index)
    a <- b + 1 
    b <- a + length(index) - 1 
    index.con[a:b] <- index
    
  } 
}

#grouping indexes in one vector index.ok
index.ok <- c(unique(index.atl[index.atl != 0]), unique(index.con[index.con != 0]))
#filtering site.locations.new with index.ok indexes 
site.locations.ok <- site.locations.new[index.ok, ]
#filtering all.species.new with site.locations.ok transects
all.species.ok <- all.species.new[transect_id %in% site.locations.ok$transect_id]