# 24 de mayo de 2022
# Aplicación piloto del flujo de trabajo: modulo geografico, modulo ambiental, modulo biogeografico

# Cargar funciones de los modulos
source("R/funciones_limpieza.R")

library(data.table)
library(sf)
library(terra)

# Explorar correlacion entre variables ambientales para ser usadas en los filtros ambientales 

envars <- list.files("envars/", ".tif", full.names = T) |> rast()
do.corr.envars(envars = envars)

#             bio_1       bio_12    radMediana   vaprCv     wc21elev_s  wc21slope
# bio_1       1.00000000  0.4048917 -0.01698513 -0.3866249 -0.97952537 -0.6453555
# bio_12      0.40489172  1.0000000 -0.21074187 -0.6297357 -0.40401289 -0.3175567
# radMediana -0.01698513 -0.2107419  1.00000000  0.1356192  0.07965861  0.1589308
# vaprCv     -0.38662488 -0.6297357  0.13561915  1.0000000  0.36276490  0.2303313
# wc21slope  -0.64535545 -0.3175567  0.15893083  0.2303313  0.67322390  1.0000000

# 1. cargar registros
load("registros/subset16.RData")
spx <- subset16[acceptedNameUsage == "Acacia decurrens"]
rm(subset16);gc()

# 2. aplicar modulo geografico: gazzeter y outliers
spx_geo <- do.geographic.label(data_base = spx, col_sp = "acceptedNameUsage", col_lon = "lon", 
                               col_lat = "lat", gazzeters = T, outliers = T)

# 3. aplicar modulo ambiental: univariado y multivariado
spx_amb <- do.environmental.label(env = envars, data_base = spx, col_lon = "lon", col_lat = "lat", 
                                  univar = T, multivar = T)

# 4. aplicar modulo biogeografico




# un monton de "" en especies name, ojo
