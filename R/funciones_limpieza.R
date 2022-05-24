# Recopilacion de funciones para desarrollar la limpieza

# A. Modulo geografico
# A.1 Etiquetado geografico

do.geographic.label <- function(data_base, col_sp,  col_lon, col_lat, gazzeters = F, outliers = F,
                                test_gazz = c("capitals", "centroids", "equal", "gbif", "institutions", 
                                              "seas", "gbif", "zeros", "duplicates"), 
                                test_ouliers = c("distance", "quantile", "mad"), 
                                outliers_details = c(1000, 5, 5), thinning_res = 0.0833, mergeto_db = F){
  
  if(!require(CoordinateCleaner)) install.packages("CoordinateCleaner")
  
  if(gazzeters == T){
    # cada test de tipo gazzeter elegido, establecido por el vector test_gazzs, se desarrolla con la 
    # función clean_coordinates() del paquete CoordianteCleaner
    gazz_results <- CoordinateCleaner::clean_coordinates(
      x = data_base,
      lon = col_lon,
      lat = col_lat,
      species = col_sp,
      tests = test_gazz,
      capitals_rad = 10000,
      centroids_rad = 1000,
      centroids_detail = "both",
      inst_rad = 100,
      range_rad = 1000,
      country_refcol = "iso_a3",
      value = "spatialvalid",
      verbose = FALSE,
      report = FALSE
    )
    gazz_results <- gazz_results[, (ncol(gazz_results)-length(test_gazz)):ncol(gazz_results)]
    gazz_results <- ((!gazz_results)*1) |> as.data.table() 
  }
  
  if(outliers == T){
    outl_list <- list()
    for(i in 1:length(test_ouliers)){
      
      # seleccionar del vector outliers details cada uno de los parametros que necesita cada test,
      # estan en orden del vector test_outliers
      if(test_ouliers[i] == "distance") tdi = outliers_details[1]
      if(test_ouliers[i] == "quantile") mltpl = outliers_details[2]
      if(test_ouliers[i] ==  "mad") mltpl = outliers_details[3]
      
      # cada test elegido, establecido por el vector test_outliers, se desarrolla con la 
      # función cc_outl() del paquete CoordianteCleaner
      outl_i <- cc_outl(
        x = data_base,
        lon = col_lon,
        lat = col_lat,
        species = col_sp,
        method = test_ouliers[i],
        tdi = tdi,
        mltpl = 5,
        value = "flagged",
        sampling_thresh = 0,
        min_occs = 7,
        thinning = T,
        thinning_res = thinning_res,
        verbose = F
      ) |> as.data.table()
      outl_list[[i]] <- outl_i
    }
    outliers_results <- do.call(cbind, outl_list)
    colnames(outliers_results) <- paste0("geo.", test_ouliers)
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

# B. Modulo ambiental

# B.1 correlaciones
do.corr.envars <- function(envars = envars, sample_size = 10000){
  if(!require(terra)) install.packages("terra")
  sample_envars <- terra::spatSample(envars, size = 10000, method = "random", replace = F, na.rm = T,
             as.df= T)
  correlacion <- sample_envars |> cor()
  return(correlacion)
}

# B.2 etiquetado ambiental
do.environmental.label <- function(data_base, col_lon, col_lat, univar = F, multivar = F,
                                   test_univar, 
                                   test_multivar, univar_details){
  
}



# generar capa espacial de puntos

gen.st.points <- function(dat, collon = col.lon, collat = col.lat) {
  st.points <- dat %>%
    dplyr::select(collon, collat) %>%
    st_as_sf(coords = c(collon, collat), crs = st_crs("EPSG:4326")) %>%
    st_transform(st_crs("EPSG:4326"))
}
