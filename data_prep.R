library(rgdal)
library(raster)
library(lubridate)
library(folderfun)
library(sf)

setff("In", "H:/Shared drives/Data/Original/pest-occurrence/Spotted Lanternfly/")

## Read in shapefile data
slf_survey_2018_19_PDA <- 
readOGR(ffIn("PA_Visual_Survey_2019_Download_1_14_2020.shp"))
slf_survey_2017 <- readOGR(ffIn("PDA_SLF_Surveys_2017_2018.shp"))
slf_2019_PPQ <- readOGR("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/raw_data/2019_AllStates_Visual_SurveyDownload_11_04_2019.shp")
slf_survery_2017_18_table <- slf_survey_2017@data
slf_survey_2019_PDA_table <- slf_survey_2018_19_PDA@data
slf_survey_2017_table <- slf_survery_2017_18_table[slf_survery_2017_18_table$ServiceYea == 2017,]
slf_survey_2018_table <- slf_survery_2017_18_table[slf_survery_2017_18_table$ServiceYea == 2018,]

names(slf_survey_2017_table)[8] <- 'Count'
names(slf_survey_2017_table)[4] <- 'Service.Date'
names(slf_survey_2017_table)[12] <- 'Primary.Surveyor'
slf_survey_2017_table$State <- 'PA'
slf_survey_2017_table$positive <- 0
slf_survey_2017_table$positive[slf_survey_2017_table$PestStatus == 'Positive'] <- 1

slf_survey_2017_ <- slf_survey_2017_table[,c("Latitude", "Longitude", "County", "State", "Count", "Host", "Primary.Surveyor", "Service.Date", "positive")]

names(slf_survey_2018_table)[8] <- 'Count'
names(slf_survey_2018_table)[4] <- 'Service.Date'
names(slf_survey_2018_table)[12] <- 'Primary.Surveyor'
slf_survey_2018_table$State <- 'PA'
slf_survey_2018_table$positive <- 0
slf_survey_2018_table$positive[slf_survey_2018_table$PestStatus == 'Positive'] <- 1

slf_survey_2018_ <- slf_survey_2018_table[,c("Latitude", "Longitude", "County", "State", "Count", "Host", "Primary.Surveyor", "Service.Date", "positive")]

slf_survey_2018sp <- readOGR("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/raw_data/SLF_survey_2018.shp")

## 2018 2019 data from shapefiles
slf_survey_2018sp2 <- spTransform(slf_survey_2018sp, CRSobj = CRS("+init=epsg:4326"))
slf_survey_2018sp2$X <- coordinates(slf_survey_2018sp2)[,1]
slf_survey_2018sp2$Y <- coordinates(slf_survey_2018sp2)[,2]
slf_survey_18_19 <- slf_survey_2018sp2@data
slf_survey_18_19$positive <- 0
slf_survey_18_19$positive[slf_survey_18_19$PestStatus == 'Positive'] <- 1
names(slf_survey_18_19)[29] <- "Longitude"
names(slf_survey_18_19)[30] <- "Latitude"
names(slf_survey_18_19)[23] <- "Primary.Surveyor"
names(slf_survey_18_19)[24] <- "Service.Date"
names(slf_survey_18_19)[11] <- "Host"
slf_survey_18_19$Year <- 2018

slf_survey_2018_ppq <- slf_survey_18_19[,c("Latitude", "Longitude", "County", "State", "Population", "Host", "Primary.Surveyor", "Service.Date", "positive", "Year")]

## 2019
slf_2019_PPQ2 <- spTransform(slf_2019_PPQ, CRSobj = CRS("+init=epsg:4326"))
slf_2019_PPQ2$X <- coordinates(slf_2019_PPQ2)[,1]
slf_2019_PPQ2$Y <- coordinates(slf_2019_PPQ2)[,2]
slf_survey_19_ppq <- slf_2019_PPQ2@data
slf_survey_19_ppq$positive <- 0
slf_survey_19_ppq$positive[slf_survey_19_ppq$PestStatus == 'Positive'] <- 1
names(slf_survey_19_ppq)[31] <- "Longitude"
names(slf_survey_19_ppq)[32] <- "Latitude"
names(slf_survey_19_ppq)[23] <- "Primary.Surveyor"
names(slf_survey_19_ppq)[29] <- "Service.Date"
names(slf_survey_19_ppq)[27] <- "Year"
names(slf_survey_19_ppq)[4] <- "Host"
slf_survey_19_ppq$Year <- 2018

slf_survey_2019_ppq <- slf_survey_19_ppq[,c("Latitude", "Longitude", "Population", "Host", "Primary.Surveyor", "Service.Date", "positive", "Year")]


#PDA
slf_survey_2018_19_PDA2 <- spTransform(slf_survey_2018_19_PDA, CRSobj = CRS("+init=epsg:4326"))
slf_survey_2018_19_PDA2$X <- coordinates(slf_survey_2018_19_PDA2)[,1]
slf_survey_2018_19_PDA2$Y <- coordinates(slf_survey_2018_19_PDA2)[,2]
slf_pda_survey_2019 <- slf_survey_2018_19_PDA2@data
slf_pda_survey_2019$positive <- 0
slf_pda_survey_2019$positive[slf_pda_survey_2019$PestStatus == "Positive"] <- 1
names(slf_pda_survey_2019)[31] <- "Longitude"
names(slf_pda_survey_2019)[32] <- "Latitude"
names(slf_pda_survey_2019)[23] <- "Primary.Surveyor"
names(slf_pda_survey_2019)[27] <- "Year"
names(slf_pda_survey_2019)[29] <- "Service.Date"
names(slf_pda_survey_2019)[4] <- "Host"
slf_pda_survey_2019$State <- 'PA'

slf_survey_2019_pda <- slf_pda_survey_2019[,c("Latitude", "Longitude", "State", "Host", "Primary.Surveyor", "Service.Date", "positive")]

## Read in Band data
slfband2015 <- read.csv(ffIn("Lyco2015_BandData_export.csv", header = TRUE)
slfband2016 <- read.csv(ffIn("Lyco2016_BandData_export.csv", header = TRUE)
slfband2017 <- read.csv(ffIn("Lyco2017_BandData_export.csv", header = TRUE)
slfband2017_treat <- read.csv(ffIn("PDA_Data_2017_Band_Treat.csv", header = TRUE)
                              
## Read in Visual Survey Data
slfvis2015 <- read.csv(ffIn("Lyco2015_VisualData_export.csv", header = TRUE)
slfvis2016 <- read.csv(ffIn("Lyco2016_VisualData_export.csv", header = TRUE)
slfvis2017 <- read.csv(ffIn("Lyco2017_VisualData_export.csv", header = TRUE)
                       
## Read in grape survey data
slfgrape2016 <- read.csv(ffIn("Lyco2016_GrapeData_export.csv", header = TRUE)
slfgrape2017 <- read.csv(ffIn("Lyco2017_GrapeData_export.csv", header = TRUE)
                         
## Read in management data
treatment_2019_new <- readOGR("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/raw_data/TreatmentPlans2019_All_States.shp")
treatment_2019 <- readOGR("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_6_state_region_psuedo_mercator/surveys.shp")

## read in plants data
plants2019 <- read.csv(ffIn("pda pa plants 2019.csv", header = TRUE)
                       
## Inaturalist data
library(spocc)
library(rinat)
species_of_interest <- "Lycorma delicatula"
extent <- c(24.396308,-124.848974, 49.384358, -66.885444) # bounding box for lower 48 states
slf_inaturalist_detections <- get_inat_obs(taxon_name = species_of_interest, maxresults = 10000, geo = TRUE, bounds = extent, quality = "research")

## Data exploration visual inspections
names(slfvis2015)
names(slfvis2016)
names(slfvis2017)

slf_2015_vis <- slfvis2015[,c("Latitude", "Longitude", "County", "State", "Count", "Order", "Family", "Genus", "Species", "Primary.Surveyor", "Service.Date", "Host", "Host.Species", "Host1", "Other.Host")]
slf_2015_vis$positive <- 0
slf_2015_vis$positive[slf_2015_vis$Count > 0 & slf_2015_vis$Genus == "Lycorma"] <- 1

slf_2016_vis <- slfvis2016[,c("Latitude", "Longitude", "County", "State", "Count", "Order", "Family", "Genus", "Species", "Primary.Surveyor", "Service.Date", "Host", "Host.Species", "Host1", "Other.Host")]
slf_2016_vis$positive <- 0
slf_2016_vis$positive[slf_2016_vis$Count > 0 & slf_2016_vis$Genus == "Lycorma"] <- 1

slf_2017_vis <- slfvis2017[,c("Latitude", "Longitude", "County", "State", "Count", "Order", "Family", "Genus", "Species", "Primary.Surveyor", "Service.Date", "Host", "Host.Species", "Host1", "Other.Host")]
slf_2017_vis$positive <- 0
slf_2017_vis$positive[slf_2017_vis$Count > 0 & slf_2017_vis$Genus == "Lycorma"] <- 1

## Data prep bands
names(slfband2015)
names(slfband2016)
names(slfband2017)
names(slfband2017_treat)

slf_2015_band <- slfband2015[,c("Latitude", "Longitude", "County", "State", "Count", "Order", "Family", "Genus", "Species", "Primary.Surveyor", "Service.Date", "Collection.Method")]
slf_2015_band$positive <- 0
slf_2015_band$positive[slf_2015_band$Count > 0 & slf_2015_band$Collection.Method == 'Sticky Band Pole'] <- 1

slf_2016_band <- slfband2016[,c("Latitude", "Longitude", "County", "State", "Count", "Order", "Family", "Genus", "Species", "Primary.Surveyor", "Service.Date", "Collection.Method")]
slf_2016_band$positive <- 0
slf_2016_band$positive[slf_2016_band$Count > 0 & slf_2016_band$Collection.Method == 'Sticky Band Pole'] <- 1

slf_2017_band <- slfband2017[,c("Latitude", "Longitude", "County", "State", "Count", "Order", "Family", "Genus", "Species", "Primary.Surveyor", "Service.Date", "Collection.Method")]
slf_2017_band$positive <- 0
slf_2017_band$positive[slf_2017_band$Count > 0 & slf_2017_band$Collection.Method == 'Sticky Band Pole'] <- 1

slf_2017_band_treat <- slfband2017_treat[,c("Latitude", "Longitude", "County", "SLF_Count", "Genus", "Species", "ServDate", "SurvTeam")]
slf_2017_band_treat$State <- "PA"
slf_2017_band_treat$positive <- 0
slf_2017_band_treat$positive[slf_2017_band_treat$SLF_Count > 0] <- 1
slf_2017_band_treat$Order <- slfband2017$Order[slfband2017$Genus == "Lycorma"][1]
slf_2017_band_treat$Family <- slfband2017$Family[slfband2017$Genus == "Lycorma"][1]
slf_2017_band_treat$Genus <- slfband2017$Family[slfband2017$Genus == "Lycorma"][1]
slf_2017_band_treat$Species <- slfband2017$Family[slfband2017$Genus == "Lycorma"][1]
slf_2017_band_treat$Count <- slf_2017_band_treat$SLF_Count
slf_2017_band_treat$Primary.Surveyor <- slf_2017_band_treat$SurvTeam
slf_2017_band_treat$Service.Date <- slf_2017_band_treat$ServDate
slf_2017_band_treat$Collection.Method <- 'Sticky Band Pole'
slf_2017_band_treat <- slf_2017_band_treat[,c("Latitude", "Longitude", "County", "State", "Count", "Order", "Family", "Genus", "Species", "Primary.Surveyor", "Service.Date", "Collection.Method", "positive")]


slf_2017_band_all <- rbind(slf_2017_band, slf_2017_band_treat)
## Data prep Grape
names(slfgrape2016)
names(slfgrape2017)

slf_2016_grape <- slfgrape2016[,c("Latitude", "Longitude", "County", "State", "Count", "Order", "Family", "Genus", "Species", "Primary.Surveyor", "Service.Date", "Habitat")]
slf_2016_grape$positive <- 0
slf_2016_grape$positive[slf_2016_grape$Count > 0 & slf_2016_grape$Genus == "Lycorma"] <- 1

slf_2017_grape <- slfgrape2017[,c("Latitude", "Longitude", "County", "State", "Count", "Order", "Family", "Genus", "Species", "Primary.Surveyor", "Service.Date", "Habitat")]
slf_2017_grape$positive <- 0
slf_2017_grape$positive[slf_2017_grape$Count > 0 & slf_2017_grape$Genus == "Lycorma"] <- 1

## SHapefile data
names(slf_survery_2017_table)

## counts from slf_detections data compared to most recent data parsed by type with negatives
counts_detections <- data.frame(table(slf_detections$Year))
names(counts_detections) <- c("year", "detections_data")
counts_detections$vis_data <- 0
counts_detections$band_data <- 0
counts_detections$grape_data <- 0

counts_detections$vis_data[counts_detections$year == 2015] <- sum(slf_2015_vis$positive)
counts_detections$vis_data[counts_detections$year == 2016] <- sum(slf_2016_vis$positive)
counts_detections$vis_data[counts_detections$year == 2017] <- sum(slf_2017_vis$positive)

counts_detections$band_data[counts_detections$year == 2015] <- sum(slf_2015_band$positive)
counts_detections$band_data[counts_detections$year == 2016] <- sum(slf_2016_band$positive)
counts_detections$band_data[counts_detections$year == 2017] <- sum(slf_2017_band_all$positive)

counts_detections$grape_data[counts_detections$year == 2016] <- sum(slf_2016_grape$positive)
counts_detections$grape_data[counts_detections$year == 2017] <- sum(slf_2017_grape$positive)

counts_detections$total <- counts_detections$vis_data + counts_detections$band_data + counts_detections$grape_data


## Bind all data sets into one
slf_2015 <- dplyr::bind_rows(slf_2015_vis, slf_2015_band)
slf_2015[slf_2015$positive > 0, c("Order", "Family", "Genus", "Species")] <- slf_2015[slf_2015$Genus == 'Lycorma', c("Order", "Family", "Genus", "Species")][1,]
slf_2015$Year <- 2015
slf_2015 <- unique(slf_2015)

slf_2016 <- dplyr::bind_rows(slf_2016_vis, slf_2016_band, slf_2016_grape)
slf_2016[slf_2016$positive > 0, c("Order", "Family", "Genus", "Species")] <- slf_2016[slf_2016$Genus == 'Lycorma', c("Order", "Family", "Genus", "Species")][1,]
slf_2016 <- slf_2016[slf_2016$Longitude < 0,]
slf_2016$Year <- 2016
slf_2016 <- unique(slf_2016)

slf_2017 <- dplyr::bind_rows(slf_2017_vis, slf_2017_band_all, slf_2017_grape, slf_survey_2017_)
slf_2017[slf_2017$positive > 0, c("Order", "Family", "Genus", "Species")] <- slf_2017[slf_2017$Genus == 'Lycorma', c("Order", "Family", "Genus", "Species")][1,]
slf_2017$Year <- 2017
slf_2017 <- unique(slf_2017)

slf_2018 <- dplyr::bind_rows(slf_survey_2018_, slf_survey_2018_ppq)
slf_2018[slf_2018$positive > 0, c("Order", "Family", "Genus", "Species")] <- slf_2017[slf_2017$Genus == 'Lycorma', c("Order", "Family", "Genus", "Species")][1,]
slf_2018$Year <- 2018
slf_2018 <- unique(slf_2018)

slf_2019 <- dplyr::bind_rows(slf_survey_2019_ppq, slf_survey_2019_pda)
slf_2019[slf_2019$positive > 0, c("Order", "Family", "Genus", "Species")] <- slf_2017[slf_2017$Genus == 'Lycorma', c("Order", "Family", "Genus", "Species")][1,]
slf_2019$Year <- 2019
slf_2019 <- unique(slf_2019)

slf_2015_2019 <- dplyr::bind_rows(slf_2015, slf_2016, slf_2017, slf_2018, slf_2019)
slf_2015_2019 <- unique(slf_2015_2019)

## Write out slf detection shapefiles 
slf_2015sp <- SpatialPointsDataFrame(slf_2015[,c(2,1)], data = slf_2015, proj4string = CRS("+init=epsg:4326"))
writeOGR(obj=slf_2015sp, dsn="H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources", layer="slf_surveys_2015", driver="ESRI Shapefile",overwrite_layer=TRUE)
writeOGR(obj=slf_2015sp, dsn="H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/slf_surveys_2015", layer="slf_surveys_2015", driver="GPKG", overwrite_layer=TRUE)

slf_2016sp <- SpatialPointsDataFrame(slf_2016[,c(2,1)], data = slf_2016, proj4string = CRS("+init=epsg:4326"))
writeOGR(obj=slf_2016sp, dsn="H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources", layer="slf_surveys_2016", driver="ESRI Shapefile",overwrite_layer=TRUE)
writeOGR(obj=slf_2016sp, dsn="H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources", layer="slf_surveys_2016", driver="ESRI Shapefile",overwrite_layer=TRUE)

slf_2017sp <- SpatialPointsDataFrame(slf_2017[,c(2,1)], data = slf_2017, proj4string = CRS("+init=epsg:4326"))
writeOGR(obj=slf_2017sp, dsn="H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources", layer="slf_surveys_2017", driver="ESRI Shapefile",overwrite_layer=TRUE)
writeOGR(obj=slf_2017sp, dsn="H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources", layer="slf_surveys_2017", driver="ESRI Shapefile",overwrite_layer=TRUE)

slf_2018sp <- SpatialPointsDataFrame(slf_2018[,c(2,1)], data = slf_2018, proj4string = CRS("+init=epsg:4326"))
writeOGR(obj=slf_2018sp, dsn="H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources", layer="slf_surveys_2018", driver="ESRI Shapefile",overwrite_layer=TRUE)
writeOGR(obj=slf_2018sp, dsn="H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/slf_surveys_2018", layer="slf_surveys_2018", driver="GPKG", overwrite_layer=TRUE)

slf_2019sp <- SpatialPointsDataFrame(slf_2019[,c(2,1)], data = slf_2019, proj4string = CRS("+init=epsg:4326"))
writeOGR(obj=slf_2019sp, dsn="H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources", layer="slf_surveys_2019", driver="ESRI Shapefile",overwrite_layer=TRUE)
writeOGR(obj=slf_2019sp, dsn="H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/slf_surveys_2019", layer="slf_surveys_2019", driver="GPKG", overwrite_layer=TRUE)

slf_2015_2019sp <- SpatialPointsDataFrame(slf_2015_2019[,c(2,1)], data = slf_2015_2019, proj4string = CRS("+init=epsg:4326"))
writeOGR(obj=slf_2015_2019sp, dsn="H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources", layer="slf_surveys_2015_2019", driver="ESRI Shapefile",overwrite_layer=TRUE)
writeOGR(obj=slf_2015_2019sp, dsn="H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/slf_surveys_2015_2019", layer="slf_surveys_2015_2019", driver="GPKG", overwrite_layer=TRUE)
