# 26 de julio de 2022
#coordenadas y nombres de archivos

library(openxlsx)
library(sf)

dataCoords <- list.files("resultadosPiloto/", ".xlsx", full.names = T)
dataNames <- list.files("resultadosPiloto/", ".xlsx") %>% gsub(".xlsx", "", x = .)

gen.st.points <- function(dat, collon = col.lon, collat = col.lat) {
  st.points <- dat %>%
    st_as_sf(coords = c(collon, collat), crs = st_crs("EPSG:4326")) %>%
    st_transform(st_crs("EPSG:4326"))
}

fol <- "resultadosPiloto/shp"
dir.create(fol, showWarnings = F)

for(i in 1:length(dataCoords)){
  a <- read.xlsx(dataCoords[i], sheet = 1)
  nm <- paste0(fol, "/", dataNames[i])
  sp_sf <- gen.st.points(a, collon = "decimalLongitude", collat = "decimalLatitude")
  write_sf(sp_sf, paste0(nm, ".gpkg"), delete_layer = T)
}


