# Recopilacion de funciones para desarrollar la limpieza

# A. Modulo geografico
# A.1 Etiquetado geografico

do.geographic.label <- function(data_base, col_sp,  col_lon, col_lat, gazzeters = F, outliers = F,
                                test_gazz = c("capitals", "centroids", "equal", "gbif", "institutions", 
                                              "seas", "gbif", "zeros", "duplicates"), 
                                test_outliers = c("dist", "iqr", "mad"), 
                                outliers_details = c("tdi" = 100, "mtpl_iqr" = 1.5, "mtpl_mad" = 5), 
                                thinning_res = 0.0833, mergeto_db = F){
  
  if(!require(CoordinateCleaner)) install.packages("CoordinateCleaner")
  
  vect <- c(col_sp, col_lon, col_lat)
  data_base <- data_base[, ..vect]
  
  if(gazzeters == T){
    # cada test de tipo gazzeter elegido, establecido por el vector test_gazzs, se desarrolla con la 
    # función clean_coordinates() del paquete CoordianteCleaner
    gazz_results <- CoordinateCleaner::clean_coordinates(
      x = data_base, lon = col_lon, lat = col_lat, species = col_sp, tests = test_gazz,
      capitals_rad = 10000, centroids_rad = 1000, centroids_detail = "both", inst_rad = 100,
      range_rad = 1000, country_refcol = "iso_a3", value = "spatialvalid", verbose = FALSE,
      report = FALSE
    )
    # generar como output unicamente las columnas de los test y eliminar el summary
    gazz_results <- gazz_results[, (ncol(gazz_results)-length(test_gazz)):(ncol(gazz_results)-1)]
    gazz_results <- ((!gazz_results)*1) |> as.data.table() 
  }
  
  if(outliers == T){
    outl_list <- list()
    
    for(i in 1:length(test_outliers)){
      
      test_i <- test_outliers[i]
      
      # seleccionar del vector outliers details cada uno de los parametros que necesita cada test
      if(test_i == "dist"){
        tdi = outliers_details["tdi"]
        method_i <- "distance"
      }else if(test_i == "iqr"){
        mltpl = outliers_details["mtpl_iqr"]
        method_i <- "quantile"    
      }else if(test_i == "mad"){
        mltpl = outliers_details["mtpl_mad"]
        method_i <- "mad"  
      }
      
      # cada test elegido, establecido por el vector test_outliers, se desarrolla con la 
      # función cc_outl() del paquete CoordianteCleaner
      outl_i <- cc_outl(
        x = data_base, lon = col_lon, lat = col_lat, species = col_sp, method = method_i,
        tdi = tdi, mltpl = mltpl, value = "flagged", sampling_thresh = 0, min_occs = 7, thinning = T,
        thinning_res = thinning_res, verbose = F) |> as.data.table()
      outl_list[[i]] <- outl_i
    }
    outliers_results <- do.call(cbind, outl_list)
    colnames(outliers_results) <- paste0("geo.", test_outliers)
    outliers_results <- ((!outliers_results)*1) |> as.data.table()
  }
  
  # compilar resultados segun la categoria de test desarollado
  if(exists("gazz_results") & !exists("outliers_results")){
     results <- gazz_results
   }else if(!exists("gazz_results") & exists("outliers_results")){
     results <- outliers_results
   }else if(exists("gazz_results") & exists("outliers_results")){
     results <- cbind(gazz_results, outliers_results)
   }
    
  return(results)
}

#-------------------------------------------------------------------------------
# B. Modulo ambiental


# B.2 etiquetado ambiental
do.environmental.label <- function(
    env = envars, data_base, col_lon, col_lat, univar = F, multivar = F,
    test_univar = c("zscore", "std", "iqr", "rjack", "mad" ), test_multivar = c("pca_error", "maha"),
    univar_details = c("thr_std" = 4, "mtpl_iqr" = 1.5, "thr_mad" = 2), 
    multivar_details = c("pval" = 0.95), min_occs = 7
    ){
  
  env_space <- gen.env.space(env. = env, data.base = data_base, col.lon = col_lon, col.lat = col_lat)
  
  if(univar == T){
    
    # inicializar output
    out_univar <- data.frame(matrix(NA, nrow = nrow(data_base), 
                                    ncol = length(test_univar)*ncol(env_space)))
    out_names <- apply(expand.grid( colnames(env_space), test_univar), 1, paste, collapse=".")
    colnames(out_univar) <- out_names
    
    if(nrow(data_base) <=  min_occs){
      return(out_univar)
      
    }else{
      
      for(i in 1:ncol(env_space)){
        
        x <-  env_space[, ..i] |> as.matrix() |> as.vector()
        xi <- na.omit(x)
        xNA <- which(is.na(x) == T)
        
        for(a in 1:length(test_univar)){
          
          test_a <- test_univar[a]
          
          if (test_a == "zscore") out_ <- zscore(x = xi )
          if (test_a == "std") out_ <- std(x = xi , threshold = univar_details["thr_std"])*1
          if (test_a == "iqr") out_ <- iqr(x = xi, mtp = univar_details["mtpl_iqr"])*1
          if (test_a == "rjack") out_ <- rjack(x = xi )
          if (test_a == "mad") out_ <- mad.(x = xi, threshold = univar_details["thr_mad"])*1
          
          name.test <- paste0(colnames(env_space)[i], ".", test_a)
          
          if(length(xNA)!= 0){
            out_ <- append(x = out_, values = rep(NA, length(xNA)), after = xNA-1 )
          }
          out_univar[, name.test] <- out_
        }
      }
    }
    
    univar_results <- out_univar |> data.table()
  }
  
  if(multivar == T){
    
    out_multivar <- data.frame(matrix(NA, nrow = nrow(data_base), ncol = length(test_multivar)))
    colnames(out_multivar) <- test_multivar
    
    # hallar los datos perdidos, no pueden ser trabajados en los test multivariados
    xNA <- apply(X = env_space, 2, FUN = function(X){which(is.na(X))}) |> unlist() |> unique()
    
    # media y escalar las variables
    mu <-  colMeans(env_space[-xNA, ])
    env_space_scale <-  scale(env_space[-xNA, ], center = mu, scale = T)
    
    if("pca_error" %in% test_multivar){
      
      #calcular PCA
      Xpca <- prcomp(x = env_space_scale)
      
      #establecer la varianza compuesta de los componentes
      summary_pca <- summary(Xpca)
      prop <- (summary_pca$sdev^2 / sum(summary_pca$sdev^2)) |> cumsum()
      
      #elegir aquellos en donde se complete el x porcentaje de varianza a trabajar

      nComp <- which(prop >= multivar_details["pval"])[1]
      
      # recrear el dataset escalado y establecer la diferencia entre el original y el predicho
      Xhat <-  Xpca$x[,1:nComp] %*% t(Xpca$rotation[,1:nComp])
      errorq_pca <- (env_space_scale - Xhat)^2
      
      # el resto de errores por cada variable estan correlacionados al 100% por lo que el
      # error de pca en una variable predicha nos habla de todas las otras
      errorq_pca <- errorq_pca[,1] |> round(5)
      
      # agregar aquellos registros con datos perdidos
      errorq_pca <- append(x = errorq_pca, values = rep(NA, length(xNA)), after = xNA-1 )
      
      out_multivar[, "pca_error"] <- errorq_pca
      
    }
    
    if("maha" %in% test_multivar){
      
      if(!require(ClassDiscovery)) install.packages("ClassDiscovery")
      
      # transponer la matriz ambiental
      env_space_scale_t <- t(env_space_scale)
      
      # generar un objeto pca sample, ya se escalaron los datos y se centraron
      spca <- SamplePCA(env_space_scale_t, usecor = F, center = F)
      
      # establecer la proporción acumulada explicada por los componentes y elegir
      # aquel en donde se da el valor de probabilidad deseado
      prop <- round(cumsum(spca@variances)/sum(spca@variances), digits=2)
      nComp <- length(prop)
      
      #MISSING: armonizar los metodos de PCA usados
      
      # establecer la distancia de mahalanobis y hallar valores de probabilidad
      maha <- mahalanobisQC(spca, nComp)
      
      # binarizar los valores de probabilidad, marcando los que esten por debajo
      # de significatividad estadistica deseada 1 - pvalor
      maha$p.value <- (maha$p.value <= (1- multivar_details["pval"]))*1 
      colnames(maha) <- c("dist.maha", "p.value")
      
      # agregar valores NA
      maha <- insertRow(existingDF = maha, newrow = rep(NA, ncol(maha)), r = xNA)
      
      # adjuntar a la data.frame de salida
      out_multivar$maha <- maha$dist.maha
      out_multivar$p.value <- maha$p.value
    }
    
    multivar_results <- out_multivar|> data.table()
  }
  
  # compilar resultados segun la categoria de test desarollado
  if(exists("univar_results") & !exists("multivar_results")){
    results <- univar_results
  }else if(!exists("univar_results") & exists("multivar_results")){
    results <- multivar_results
  }else if(exists("univar_results") & exists("multivar_results")){
    results <- cbind(univar_results, multivar_results)
  }
  
  return(results)

}

#-------------------------------------------------------------------------------
# C.1 Modulo biogeografico
do.biogeographic.label <- function(
    data_base, col_sp, col_lon, col_lat, col_dup, shapeLayers = NULL, polynames,
    test_biog = c("bracatus", "speciesGeo", "IaVH"), biog_details = c("pval" = 0.95),
    min_occs = 7
){
  
  out_biogeo <- list()
  
  sf::sf_use_s2(FALSE)
  
  data_base <- as.data.frame(data_base)
  
  #Extraer columnas de especie y coordenadas de la base de datos de registros de presencia
  pts <- data_base[, c(col_sp, col_lon, col_lat)]
  rm(data_base);gc()
  
  # Loop que clasifica cada capa shape que fue cargado
  for(i in 1:length(shapeLayers)){
    print(paste0("shapelayer = ", i))
    if("speciesGeo" %in% test_biog){    
      
      shp_nm <- abbreviate(names(shapeLayers)[i])
      nm <- paste0("SpGeo", ".", shp_nm)
      
      #extraer de la lista de shapefiles el indice i
      shapeLayers_i <- shapeLayers[[i]]
      polynames_i <- polynames[i]
      
      # darle nombres al data frame de la base de datos aceptados por species geocoder
      colnames(pts) <- c("species", "decimalLongitude", "decimalLatitude")
      
      # extraer poligonos en donde caen los registros de presencia
      shp_ref <- shapeLayers_i[gen.st.points(pts, collon = "decimalLongitude",
                                             collat = "decimalLatitude"), ]
      

      
      # clasificar los puntos en los poligonos del shapefile de referencia usando una frecuencia de
      # 1-pval en terminos de porcentaje
      sp.class <- SpGeoCod(x = pts, y = as_Spatial(shp_ref), areanames = polynames_i, 
                           occ.thresh =  (1-biog_details["pval"])*100)
      
      # extraer aquellos poligonos en donde quedaron clasificados los registros
      geocode <- sp.class$spec_table
      geocode <- colnames(geocode)[which(geocode[1, ] != 0)]
      
      # Extraer la clasificación de cada registro en el poligono
      samples <- sp.class$samples[, "homepolygon"]
      
      # Transformar la clasificación: aquellos registros clasificados darle valor de 0,
      # aquellos que no fueron clasificados dado que estan dentro del margen de error
      # darle valor 1
      geocode_samples <- rep(NA, length(samples))
      for(x in 1:length(geocode)){
        geocode_x <- which(samples == geocode[x])
        geocode_samples[geocode_x] <- 0
      }
      geocode_samples[is.na(geocode_samples)] <- 1
      
      # Crear tabla de resultados de speciesGeocodeR
      results_SpeGeo <- data.frame(geocode_samples)
      colnames(results_SpeGeo) <- nm
      
      # agregar al shapefile una columna que identifique cada poligono como anomalo o no
      shp_polyname_vect <- st_drop_geometry(shp_ref) %>% dplyr::pull(polynames_i)
      shp_ref[, nm] <- shp_polyname_vect %in% geocode
    }  
    
    if("bracatus" %in% test_biog){
      if(!("speciesGeo" %in% test_biog) & !("IaVH" %in% test_biog)){
        stop("Si desea correr bracatus es necesario que esten activos alguno de los dos metodos de 
           clasificación de puntos inicial: 'speciesGeo' o 'IaVH'")
      }
      shp_ref_sp <- shp_ref %>% as_Spatial()
      logic_ref <- shp_ref_sp@data[ , nm] %>% unique()
      
      if(length(logic_ref) >= 2){
        range_map <- rangeMaps(range = shp_ref_sp, biogeo = nm, native = T, alien = F)
        signals <- signalCalculation (ref_reg = range_map, pts = pts, biogeo = F)
        acc <- accuracy (signals) %>% dplyr::select(accuracy)
      }else{
        acc <- data.frame(rep(1, nrow(pts)))
      }
      results_bracatus <- acc
      colnames(results_bracatus) <- paste0("bRa", ".", shp_nm)  
    }
    
    # compilar resultados segun la categoria de test desarollado
    if(exists("results_SpeGeo") & !exists("results_bracatus")){
      results_i <- results_SpeGeo
    }else if(exists("results_SpeGeo") & exists("results_bracatus")){
      results_i <- cbind(results_SpeGeo, results_bracatus)
    }
    
    out_biogeo[[i]] <- results_i
  }

  return(dplyr::bind_cols(out_biogeo))
  
}


#-------------------------------------------------------------------------------
# Funciones individuales
#-------------------------------------------------------------------------------
# Correlaciones ambientales
do.corr.envars <- function(envars = envars, sample_size = 10000){
  if(!require(terra)) install.packages("terra")
  sample_envars <- terra::spatSample(envars, size = 10000, method = "random", replace = F, na.rm = T,
                                     as.df= T)
  correlacion <- sample_envars |> cor()
  return(correlacion)
}

#-------------------------------------------------------------------------------
# Extraer espacio ambiental
gen.env.space <- function(env. = env, data.base = data_base, col.lon = col_lon, col.lat = col_lat){
  vect <- c(col.lon, col.lat)
  geodata <- data.base[, ..vect] |> data.frame() 
  env.space <- terra::extract(env., geodata)[, -1] |> data.table()
  return(env.space)
}

#-------------------------------------------------------------------------------
# generar capa espacial de puntos
gen.st.points <- function(dat, collon = col.lon, collat = col.lat) {
  st.points <- dat |>
    dplyr::select(collon, collat) |>
    st_as_sf(coords = c(collon, collat), crs = st_crs("EPSG:4326")) |>
    st_transform(st_crs("EPSG:4326"))
}
#-------------------------------------------------------------------------------
# zscore: Regi0

zscore <- function(x){
  values <- (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
  return( values )
}

#-------------------------------------------------------------------------------
# desviacion estandar: Regi0

std <- function(x, threshold){
  std <- sd(x, na.rm = T)
  mean <- mean(x, na.rm = T)
  return ( (x < mean - (threshold * std)) | (x > mean + (threshold * std)) )
}

#-------------------------------------------------------------------------------
# desviacion mediana absoluta

mad. <- function(x, threshold){
  madx <- stats::mad(x, na.rm = T)
  medianx <- median(x, na.rm = T)
  return ( (x < medianx - (threshold * madx)) | (x > medianx + (threshold * madx)) )
}
#-------------------------------------------------------------------------------
# rango intercuartilico: Regi0

iqr <- function(x, mtp){
  rangeIQR <- IQR(x, na.rm = T)
  quants <- quantile(x, na.rm = T)
  q1 <- quants[2]
  q3 <- quants[4]
  return ( (x < q1 - (mtp * rangeIQR)) | (x > q3 + (mtp * rangeIQR)) ) 
}

#-------------------------------------------------------------------------------
# rjackknife: biogeo

rjack <- function (x) 
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

insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}
