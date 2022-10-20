# C.1 Modulo biogeografico

do_biogeographic_tag <- function(
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
      shp_ref <- shapeLayers_i[gen_sf_points(pts, collon = "decimalLongitude",
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
        geocode_samples[geocode_x] <- 1
      }
      geocode_samples[is.na(geocode_samples)] <- 0
      
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
      
      if(!require(raster)) install.packages("raster")
      
      shp_ref_sp <- shp_ref %>% as_Spatial()
      logic_ref <- shp_ref_sp@data[ , nm] %>% unique()
      
      if(length(logic_ref) >= 2){
        range_map <- range_maps(range = shp_ref_sp, biogeo = nm, native = T, alien = F)
        signals <- signalCalculation (ref_reg = range_map, pts = pts, biogeo = F) #time consuming
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