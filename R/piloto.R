# 24 de mayo de 2022
# Aplicación piloto del flujo de trabajo: modulo geografico, modulo ambiental, modulo biogeografico

# Cargar funciones de los modulos
source("R/funciones_limpieza.R")

library(data.table)
library(sf)
library(terra)
library(speciesgeocodeR)
library(bRacatus)
library(dplyr)
library(openxlsx)

#--------------------
# Explorar correlacion entre variables ambientales para ser usadas en los filtros ambientales 

envars <- list.files("envars/", ".tif", full.names = T) |> rast()
do.corr.envars(envars = envars)

#             bio_1       bio_12    radMediana   vaprCv     wc21elev_s  wc21slope
# bio_1       1.00000000  0.4048917 -0.01698513 -0.3866249 -0.97952537 -0.6453555
# bio_12      0.40489172  1.0000000 -0.21074187 -0.6297357 -0.40401289 -0.3175567
# radMediana -0.01698513 -0.2107419  1.00000000  0.1356192  0.07965861  0.1589308
# vaprCv     -0.38662488 -0.6297357  0.13561915  1.0000000  0.36276490  0.2303313
# wc21slope  -0.64535545 -0.3175567  0.15893083  0.2303313  0.67322390  1.0000000

#--------------------
# leer las capas biogeograficas que se desean trabajar
hybas3 <- sf::read_sf("biogeographic/hybas_sa_lev03_v1c.shp") 
morrone <- sf::read_sf("biogeographic/Lowenberg_Neto_2014_biomod.shp") 
ecoreg <-  sf::read_sf("biogeographic/wwf_terr_ecos_biomod.shp") 
bior <- sf::read_sf("biogeographic/bior164.shp")

# generar una lista y darle nombre a cada capa en la lista
list_bio <- list(hybas3, morrone, ecoreg, bior); names(list_bio) <- c("hybas3", "morrone", "ecoreg", "bior")

#--------------------
# 1. cargar registros
spx <- fread("registros/RegistrosBM19042022.csv", stringsAsFactors = F) |> 
  filter(spatialDuplicated == "false" ) %>% mutate("decimalLongitude" = as.numeric(decimalLongitude),
                                                   "decimalLatitude" = as.numeric(decimalLatitude))

# 2.etiqueta de experto
experto <- (is.na(spx$reportedDate))*1
tabexp <- table(spx$acceptedNameUsage, experto) %>% data.frame()
write.xlsx(tabexp, "resultadosPiloto2.xlsx",sheet = 1)

spxList <- split(spx,f = spx$acceptedNameUsage)

# Elegir especies que al menos tengan 10 registros y con mas de 1 curado por experto
for(i in 1:length(spxList)){
  if((tabexp[tabexp$experto == 0, 3][i]) == 0 | ((tabexp[tabexp$experto == 1, 3][i]) <= 10)){
    spxList[[i]] <- NA    
  }
}
spxList <- spxList[-which(is.na(spxList))]

# carpeta para escribir tablas
fol <- "resultadosPiloto"
dir.create(fol, showWarnings = F)

#--------------------
for(i in 251:length(spxList)){
  spx_i <- spxList[[i]]
    spp_i <- spx_i$acceptedNameUsage[1]
    print(spp_i)
    
    filelog <- file(paste0(fol,"/", gsub(" ", "_", spp_i), ".txt"), "w")
    
    time1 <-Sys.time()
    
    # aplicar modulo geografico: gazzeter y outliers
    spx_geo <- do.geographic.label(data_base = spx_i, col_sp = "acceptedNameUsage", col_lon = "decimalLongitude", 
                                   col_lat = "decimalLatitude", gazzeters = T, outliers = T)
    
    time2 <- give.msg.time(time.1 = time1)
    
    writeLines(paste0("Geografico: ", "\n", time2), filelog)
    
    # aplicar modulo ambiental: univariado y multivariado
    spx_amb <- do.environmental.label(env = envars, data_base = spx_i, col_lon = "decimalLongitude", 
                                      col_lat = "decimalLatitude", univar = T, multivar = T)
    
    time3 <- give.msg.time(time.1 = time1)
    writeLines(paste0("Ambiental: ", "\n", time3), filelog)
    
    # aplicar modulo biogeografico
    
    spx_biog <- do.biogeographic.label(data_base = spx_i, col_sp = "acceptedNameUsage", col_lon = "decimalLongitude", 
                                       col_lat = "decimalLatitude", shapeLayers = list_bio, 
                                       polynames = c("HYBAS_ID", "NUMB", "BIOME", "BIOME_NAME"),
                                       test_biog = c("speciesGeo"),
                                       biog_details = c("pval" = 0.95))
    
    time4 <- give.msg.time(time.1 = time1)
    writeLines(paste0("Ambiental: ", "\n", time4), filelog)
    
    close(filelog)
    # tabla general
    a <- data.frame("acceptedNameUsage" = spx_i$acceptedNameUsage, "decimalLatitude" = spx_i$decimalLatitude, 
                    "decimalLongitude" = spx_i$decimalLongitude, spx_geo, spx_amb, spx_biog, 
                    "experto" = is.na(spx_i$reportedDate)*1)
    write.xlsx(a, paste0("resultadosPiloto/", spp_i,".xlsx"), rowNames = F)
}

# Marinas: datos fuera de las variables usadas
# "Pristimantis merostictus"
# Error in solve.default(covMat) : 
#   system is computationally singular: reciprocal condition number = 1.06734e-19 

# [1] "Rhinoclemmys annulata"
# Error in CoordinateCleaner::clean_coordinates(x = data_base, lon = col_lon,  : 
# invalid coordinates found in rows, clean dataset before proceeding:
# 20 
# 21 
# 28 
# 34 
# 38 
# In addition: There were 50 or more warnings (use warnings() to see the first 50)

# [1] "Thryophilus sernai"
# Error in CoordinateCleaner::clean_coordinates(x = data_base, lon = col_lon,  : 
# invalid coordinates found in rows, clean dataset before proceeding:
# 80 
# 81 
# 84 
# 106 
# 129 
# 133 
# 136 
# 138 
# In addition: There were 50 or more warnings (use warnings() to see the first 50)