#-------------------------------------------------------------------------------
# Funciones individuales
#-------------------------------------------------------------------------------
# Correlaciones ambientales
do_corr_envars <- function(envars = envars, sample_size = 10000){
  if(!require(terra)) install.packages("terra")
  sample_envars <- terra::spatSample(envars, size = 10000, method = "random", replace = F, na.rm = T,
                                     as.df= T)
  correlacion <- sample_envars |> cor()
  return(correlacion)
}

#-------------------------------------------------------------------------------
# Extraer espacio ambiental
gen_env_space <- function(env. = env, data.base = data_base, col.lon = col_lon, col.lat = col_lat){
  vect <- c(col.lon, col.lat)
  geodata <- data.base[, ..vect] |> data.frame() 
  env.space <- terra::extract(env., geodata)[, -1] |> data.table()
  return(env.space)
}

#-------------------------------------------------------------------------------
# generar capa espacial de puntos
gen_sf_points <- function(dat, collon = col.lon, collat = col.lat) {
  st.points <- dat |>
    dplyr::select(collon, collat) |>
    st_as_sf(coords = c(collon, collat), crs = st_crs("EPSG:4326")) |>
    st_transform(st_crs("EPSG:4326"))
}
#-------------------------------------------------------------------------------
# zscore: Regi0

calc_zscore <- function(x){
  values <- (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
  return( values )
}

#-------------------------------------------------------------------------------
# desviacion estandar: Regi0

calc_std <- function(x, threshold){
  std <- sd(x, na.rm = T)
  mean <- mean(x, na.rm = T)
  return ( (x < mean - (threshold * std)) | (x > mean + (threshold * std)) )
}

#-------------------------------------------------------------------------------
# desviacion mediana absoluta

calc_mad <- function(x, threshold){
  madx <- stats::mad(x, na.rm = T)
  medianx <- median(x, na.rm = T)
  return ( (x < medianx - (threshold * madx)) | (x > medianx + (threshold * madx)) )
}
#-------------------------------------------------------------------------------
# rango intercuartilico: Regi0

calc_iqr <- function(x, mtp){
  rangeIQR <- IQR(x, na.rm = T)
  quants <- quantile(x, na.rm = T)
  q1 <- quants[2]
  q3 <- quants[4]
  return ( (x < q1 - (mtp * rangeIQR)) | (x > q3 + (mtp * rangeIQR)) ) 
}

#-------------------------------------------------------------------------------
# rjackknife: biogeo

calc_rjack <- function (x) 
{
  xout <- x
  
  xx <- x
  x <- unique(x)
  rng <- diff(range(x, na.rm = T))
  mx <- mean(x, na.rm = T)
  n <- as.numeric(length(x))
  n1 <- abs(n - 1)
  t1 <- (0.95 * sqrt(n)) + 0.2
  x <- sort(x)
  y <- rep(0, n1)
  
  for (i in 1:n1) {
    x1 <- x[i + 1]
    if (x[i] < mx) {
      y[i] <- (x1 - x[i]) * (mx - x[i])
    }else {
      y[i] <- (x1 - x[i]) * (x1 - mx)
    }
  }
  
  my <- mean(y, na.rm = T)
  z <- y/(sqrt(sum((y - my)^2)/n1))
  out <- rep(0, length(xx))
  
  if (any(z > t1)) {
    f <- which(z > t1)
    v <- x[f]
    if (v < median(x)) {
      xa <- (xx <= v) * 1
      out <- out + xa
    }
    if (v > median(x)) {
      xb <- (xx >= v) * 1
      out <- out + xb
    }
  }else{
    out <- out
  }
  
  return(out)
}


#-------------------------------------------------------------------------------
# insertar fila en una posición especifica sin dependen de paquetes adicionales
# https://stackoverflow.com/questions/11561856/add-new-row-to-dataframe-at-specific-row-index-not-appended

insert_row <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

#-------------------------------------------------------------------------------
range_maps <- function (range, biogeo = "legend", native = "Extant (resident)", 
                        alien = "Introduced") 
{
  range_native <- range[which(range@data[, biogeo] %in% native), 
  ]
  range_alien <- range[which(range@data[, biogeo] %in% alien), 
  ]
  
  raster_2degrees <- raster::raster(vals = NA, res = 2)
  raster_cut <- raster::crop(raster_2degrees, raster::extent(range) + c(-2, 2, -2, 2))
  range2 <- rasterize(range, raster_cut, getCover = TRUE)
  range3 <- raster::rasterToPolygons(range2)
  range4 <- range3[which(range3$layer != 0), ]
  range4$data <- 1
  range4 <- range4[which(range4$layer >= 0.05), ]
  range4$area <- raster::area(range4)/1e+06
  range_native2 <- rasterize(range_native, raster_cut, getCover = TRUE)
  range_native3 <- raster::rasterToPolygons(range_native2)
  range_native4 <- range_native3[which(range_native3$layer != 
                                         0), ]
  range_native4$data <- 1
  range_native4 <- range_native4[which(range_native4$layer >= 
                                         0.05), ]
  range_native4@data
  range_native4$area <- raster::area(range_native4)/1e+06
  range_alien2 <- rasterize(range_alien, raster_cut, getCover = TRUE)
  range_alien3 <- raster::rasterToPolygons(range_alien2)
  range_alien4 <- range_alien3[which(range_alien3$layer != 
                                       0), ]
  range_alien4$data <- 1
  range_alien4 <- range_alien4[which(range_alien4$layer >= 0.01), ]
  
  range_alien4$area <- raster::area(range_alien4)/1e+06
  range_list <- list(range4, range_native4, range_alien4)
  names(range_list) <- c("Presence", "Native", 
                         "Alien")
  return(range_list)
}
#--------------------
give_msg_time <- function(time.1) {
  time.2 <- Sys.time()
  
  # Time spent
  
  alltime <- round(as.numeric(difftime(time.2, time.1, units = "mins")), 2)
  
  t.sec <- round((alltime - trunc(alltime)) * 60)
  
  msg <- paste("Completed in ", trunc(alltime), "min", t.sec, "sec\n")
  
  return(msg)
}
