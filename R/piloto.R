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

#               bio_1      bio_12    radMediana  tminCv     vaprCv     vaprMediana  wc21elev_s  wc21slope
# bio_1        1.00000000  0.4202742 -0.01800292 -0.7073483 -0.4168254   0.9567848 -0.97955850 -0.6643378
# bio_12       0.42027415  1.0000000 -0.22197282 -0.3761027 -0.6497366   0.5257991 -0.41365203 -0.3368516
# radMediana  -0.01800292 -0.2219728  1.00000000  0.0107829  0.1333839  -0.1444367  0.07954542  0.1430822
# wc21elev_s  -0.97955850 -0.4136520  0.07954542  0.6990616  0.3924701  -0.9601748  1.00000000  0.6917267

# Aunque la elevación tiene una gran correlacion se mantendra porque para muchas especies se tienen
# acotadas las zonas de elevaciones en las que viven, asi que se mantiene

# 1. cargar registros
load("registros_demo/subset16.RData")
spx <- subset16[acceptedNameUsage == "Acacia decurrens"]
rm(subset16);gc()

# 2. aplicar modulo geografico: gazzeter y outliers
spx_geo <- do.geographic.label(data_base = spx, col_sp = "acceptedNameUsage", col_lon = "lon", col_lat = "lat",
                         gazzeters = F, outliers = T)

# 3. aplicar modulo ambiental: univariado y multivariado
spx_amb <- do.environmental.label()










# un monton de "" en especies name, ojo
