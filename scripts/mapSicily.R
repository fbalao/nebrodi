
library(elevatr)
library(rgeoboundaries)
library(sf)
library(raster)
library(ggplot2)
library(viridis)
library(rnaturalearth)
library(rnaturalearthhires)

swiss_bound <- rgeoboundaries::geoboundaries("Italy")
elevation_data <- elevatr::get_elev_raster(locations = swiss_bound, z = 9, clip = "locations")


elevation_data <- as.data.frame(elevation_data, xy = TRUE)
colnames(elevation_data)[3] <- "elevation"
# remove rows of data frame with one or more NA's,using complete.cases
elevation_data <- elevation_data[complete.cases(elevation_data), ]

ggplot() +
  geom_raster(data = elevation_data, aes(x = x, y = y, fill = elevation)) +
  geom_sf(data = swiss_bound, color = "white", fill = NA) +
  coord_sf(xlim = c(11.08, 15.13), ylim = c(36.26, 38.53)) +
  scale_fill_viridis_c() +
  labs(title = "", x = "Longitude", y = "Latitude", fill = "Elevation (meters)")




library(rnaturalearth)
library(rnaturalearthhires)

sf_swiss <- ne_states(country = "italy", returnclass = "sf")

elevation_1 <- elevatr::get_elev_raster(locations = sf_swiss, z = 7, clip = "locations")
cropped_elev <- crop(elevation_1, sf_swiss)
elevate <- as.data.frame(cropped_elev, xy = TRUE)

colnames(elevate)[3] <- "elevation_value"
elevate <- elevate[complete.cases(elevate), ]
elevate [elevate$elevation_value<0,3] <-0
library(tidyverse)
elevate <- elevate %>% mutate(across(c(x, y), round, digits = 6 ))

ggplot() +
  geom_sf(data = st_as_sfc(st_bbox(elevation_1)),color = "grey", fill = "grey",alpha = 0.05) +
  geom_raster(data = elevate, aes(x = x, y = y, fill = elevation_value)) +
  #geom_sf(data = sf_swiss, color = "grey50", fill = NA) +
  coord_sf(xlim = c(12.1, 15.6), ylim = c(36.6, 38.33)) +
  scale_fill_gradientn(colours = grey.colors(50), name = "elevation") +
  #scale_fill_viridis_c() +
  labs(title = "", x = "Longitude", y = "Latitude", fill = "Elevation (meters)")

