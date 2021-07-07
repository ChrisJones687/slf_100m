library(rgdal)
library(lubridate)
library(folderfun)
library(sf)
library(terra)
## weather data set up
setff("weather", "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/whole_usa/")
setff("out", "H:/Shared drives/Data/Raster/Regional/SLF_100m/")
slf_toh <- rast(ffout("tree_of_heaven_100m.tif"))

ct_2015 <- rast(ffweather("temperature_crit_slf_usa_2015.tif"))
ct_2015 <- project(ct_2015, slf_toh)
ct_2015 <- crop(ct_2015, slf_toh)
ct_2015 <- resample(ct_2015, slf_toh)

ct_2016 <- rast(ffweather("temperature_crit_slf_usa_2016.tif"))
ct_2016 <- project(ct_2016, slf_toh)
ct_2016 <- crop(ct_2016, slf_toh)
ct_2016 <- resample(ct_2016, slf_toh)

ct_2017 <- rast(ffweather("temperature_crit_slf_usa_2017.tif"))
ct_2017 <- project(ct_2017, slf_toh)
ct_2017 <- crop(ct_2017, slf_toh)
ct_2017 <- resample(ct_2017, slf_toh)

ct_2018 <- rast(ffweather("temperature_crit_slf_usa_2018.tif"))
ct_2018 <- project(ct_2018, slf_toh)
ct_2018 <- crop(ct_2018, slf_toh)
ct_2018 <- resample(ct_2018, slf_toh)

ct_2019 <- rast(ffweather("temperature_crit_slf_usa_2019.tif"))
ct_2019 <- project(ct_2019, slf_toh)
ct_2019 <- crop(ct_2019, slf_toh)
ct_2019 <- resample(ct_2019, slf_toh)

ct_2020 <- ct_2019

writeRaster(ct_2015, ffout("crit_temp_2015.tif"), overwrite = TRUE)
writeRaster(ct_2016, ffout("crit_temp_2016.tif"), overwrite = TRUE)
writeRaster(ct_2017, ffout("crit_temp_2017.tif"), overwrite = TRUE)
writeRaster(ct_2018, ffout("crit_temp_2018.tif"), overwrite = TRUE)
writeRaster(ct_2019, ffout("crit_temp_2019.tif"), overwrite = TRUE)
writeRaster(ct_2020, ffout("crit_temp_2020.tif"), overwrite = TRUE)

## old weather data
temp_2015 <- rast(ffweather("temp_coeff_slf_usa_2015.tif"))
temp_2015 <- project(temp_2015, slf_toh)
temp_2015 <- crop(temp_2015, slf_toh)
temp_2015 <- resample(temp_2015, slf_toh)

temp_2016 <- rast(ffweather("temp_coeff_slf_usa_2016.tif"))
temp_2016 <- project(temp_2016, slf_toh)
temp_2016 <- crop(temp_2016, slf_toh)
temp_2016 <- resample(temp_2016, slf_toh)

temp_2017 <- rast(ffweather("temp_coeff_slf_usa_2017.tif"))
temp_2017 <- project(temp_2017, slf_toh)
temp_2017 <- crop(temp_2017, slf_toh)
temp_2017 <- resample(temp_2017, slf_toh)

temp_2018 <- rast(ffweather("temp_coeff_slf_usa_2018.tif"))
temp_2018 <- project(temp_2018, slf_toh)
temp_2018 <- crop(temp_2018, slf_toh)
temp_2018 <- resample(temp_2018, slf_toh)

temp_2019 <- rast(ffweather("temp_coeff_slf_usa_2019.tif"))
temp_2019 <- project(temp_2019, slf_toh)
temp_2019 <- crop(temp_2019, slf_toh)
temp_2019 <- resample(temp_2019, slf_toh)

temp_2020 <- temp_2019

writeRaster(temp_2015, ffout("temp_2015.tif"), overwrite = TRUE)
writeRaster(temp_2016, ffout("temp_2016.tif"), overwrite = TRUE)
writeRaster(temp_2017, ffout("temp_2017.tif"), overwrite = TRUE)
writeRaster(temp_2018, ffout("temp_2018.tif"), overwrite = TRUE)
writeRaster(temp_2019, ffout("temp_2019.tif"), overwrite = TRUE)
writeRaster(temp_2020, ffout("temp_2020.tif"), overwrite = TRUE)

temp <- c(temp_2017, temp_2018, temp_2019)
crit <- c(ct_2017, ct_2018, ct_2019)

writeRaster(crit, ffout("crit_temp.tif"), overwrite = TRUE)
writeRaster(temp, ffout("temp.tif"), overwrite = TRUE)


## new weather data
setff("weather", "H:/My Drive/Folder/")
temp_2015 <- rast(ffweather("slf_temp_coeff_2015.tif"))
temp_2015 <- project(temp_2015, slf_toh)
temp_2015 <- crop(temp_2015, slf_toh)
temp_2015 <- resample(temp_2015, slf_toh)

temp_2016 <- rast(ffweather("slf_temp_coeff_2016.tif"))
temp_2016 <- project(temp_2016, slf_toh)
temp_2016 <- crop(temp_2016, slf_toh)
temp_2016 <- resample(temp_2016, slf_toh)

temp_2017 <- rast(ffweather("slf_temp_coeff_2017.tif"))
temp_2017 <- project(temp_2017, slf_toh)
temp_2017 <- crop(temp_2017, slf_toh)
temp_2017 <- resample(temp_2017, slf_toh)

temp_2018 <- rast(ffweather("slf_temp_coeff_2018.tif"))
temp_2018 <- project(temp_2018, slf_toh)
temp_2018 <- crop(temp_2018, slf_toh)
temp_2018 <- resample(temp_2018, slf_toh)

temp_2019 <- rast(ffweather("slf_temp_coeff_2019.tif"))
temp_2019 <- project(temp_2019, slf_toh)
temp_2019 <- crop(temp_2019, slf_toh)
temp_2019 <- resample(temp_2019, slf_toh)

temp_2020 <- temp_2019

writeRaster(temp_2015, ffout("new_weather/temp_2015.tif"), overwrite = TRUE)
writeRaster(temp_2016, ffout("new_weather/temp_2016.tif"), overwrite = TRUE)
writeRaster(temp_2017, ffout("new_weather/temp_2017.tif"), overwrite = TRUE)
writeRaster(temp_2018, ffout("new_weather/temp_2018.tif"), overwrite = TRUE)
writeRaster(temp_2019, ffout("new_weather/temp_2019.tif"), overwrite = TRUE)
writeRaster(temp_2020, ffout("new_weather/temp_2020.tif"), overwrite = TRUE)

temp <- c(temp_2017, temp_2018, temp_2019)

writeRaster(temp, ffout("new_weather/temp.tif"), overwrite = TRUE)

setff("In", "H:/Shared drives/Data/Raster/Regional/SLF_1km/new_weather/")
t1 <- rast(ffIn("temp_2020.tif"))
t2 <- rast(ffIn("temp_2015.tif"))
t3 <- rast(ffIn("temp_2016.tif"))
t4 <- rast(ffIn("temp_2019.tif"))
t5 <- rast(ffIn("temp_2017.tif"))
temp5 <- c(t1, t2, t3, t4, t5)
terra::writeRaster(temp5, ffIn("temp5years.tif"))
 
c1 <- rast(ffIn("temp_2020.tif"))
c2 <- rast(ffIn("temp_2015.tif"))
c3 <- rast(ffIn("temp_2016.tif"))
c4 <- rast(ffIn("temp_2019.tif"))
c5 <- rast(ffIn("temp_2017.tif"))
ctemp5 <- c(c1, c2, c3, c4, c5)

terra::writeRaster(ctemp5, ffIn("crit_temp5years.tif"))


