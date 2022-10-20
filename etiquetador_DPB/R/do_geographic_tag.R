# A. Modulo geografico
# A.1 Etiquetado geografico

do_geographic_tag <- function(data_base, col_sp,  col_lon, col_lat, gazzeters = F, outliers = F,
                              test_gazz = c("capitals", "centroids", "equal", "gbif", "institutions", 
                                            "seas", "gbif", "zeros"), 
                              test_outliers = c("dist", "iqr", "mad"), 
                              outliers_details = list("tdi" = c(300,400, "variog"), "mtpl_iqr" = c(1.5, 3, 5), 
                                                      "mtpl_mad" = c(1.5, 3, 5)), 
                              thinning_res = 0.0833, mergeto_db = F, rast_DEM ="envars/wc21elev_s.tif"){
  
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
    gazz_results <- ((gazz_results)*1) |> as.data.table() 
  }
  
  if(outliers == T){
    outl_list <- list()
    
    
    for(i in 1:length(test_outliers)){
      
      test_i <- test_outliers[i]
      
      # seleccionar de la lista outliers_details cada uno de los parametros que necesita cada test
      if(test_i == "dist"){
        tdi = outliers_details[["tdi"]]
        method_i <- "distance"
      }else if(test_i == "iqr"){
        mltpl = outliers_details[["mtpl_iqr"]]
        method_i <- "quantile"    
      }else if(test_i == "mad"){
        mltpl = outliers_details[["mtpl_mad"]]
        method_i <- "mad"  
      }
      
      if(test_i == "iqr"){
        outl_mltpl <- list()
        for(a in 1:length(mltpl)){
          outl_a <- cc_outl(
            x = data_base, lon = col_lon, lat = col_lat, species = col_sp, method = method_i,
            tdi = tdi, mltpl = mltpl[a], value = "flagged", sampling_thresh = 0, min_occs = 7, thinning = T,
            thinning_res = thinning_res, verbose = F) |> as.data.table()  
          outl_mltpl[[a]] <- outl_a
        }
        outl_mltpl <- do_call(cbind, outl_mltpl)
        colnames(outl_mltpl) <- paste0("_", mltpl)
        
        outl_list[["iqr"]] <- outl_mltpl
      }
      
      if(test_i == "mad"){
        outl_mltpl <- list()
        for(a in 1:length(mltpl)){
          outl_a <- cc_outl(
            x = data_base, lon = col_lon, lat = col_lat, species = col_sp, method = method_i,
            tdi = tdi, mltpl = mltpl[a], value = "flagged", sampling_thresh = 0, min_occs = 7, thinning = T,
            thinning_res = thinning_res, verbose = F) |> as.data.table()  
          outl_mltpl[[a]] <- outl_a
        }
        outl_mltpl <- do_call(cbind, outl_mltpl)
        colnames(outl_mltpl) <- paste0("_", mltpl)
        
        outl_list[["mad"]] <- outl_mltpl
      }
      
      if(test_i == "dist"){
        index_num <- grepl("(\\d+)", x = tdi)
        if(sum(index_num) != 0){
          num_tdi <- tdi[index_num] |> as.numeric()
          outl_tdi <- list()
          for(a in 1:length(num_tdi)){
            outl_a <- cc_outl(
              x = data_base, lon = col_lon, lat = col_lat, species = col_sp, method = method_i,
              tdi = num_tdi[a], mltpl = mltpl, value = "flagged", sampling_thresh = 0, min_occs = 7, thinning = T,
              thinning_res = thinning_res, verbose = F) |> as.data.table()
            outl_tdi[[a]] <- outl_a
          }
          outl_tdi <- do_call(cbind, outl_tdi)
          colnames(outl_tdi) <- paste0("_", num_tdi, "_km")
          
          outl_list[["tdi"]] <- outl_tdi
          
        }
        
        if("variog" %in% tdi){
          if(!require(automap)) install.packages("automap")
          if(!require(sp)) install.packages("sp")
          
          r <- rast_DEM %>% rast()
          tmp <- data_base
          tmp$r_DEM <- extract(r, as.data.frame(data_base)[,c(col_lon, col_lat)])[,2]
          tmp <- na.omit(tmp)
          coordinates(tmp) <-  ~ decimalLongitude + decimalLatitude
          emp_Var <- autofitVariogram(formula = r_DEM ~ 1, input_data =  tmp)
          
          distx <- emp_Var$var_model[2,"range"]*111
          
          outl_variog <- cc_outl(
            x = data_base, lon = col_lon, lat = col_lat, species = col_sp, method = method_i,
            tdi = distx, mltpl = mltpl, value = "flagged", sampling_thresh = 0, min_occs = 7, thinning = T,
            thinning_res = thinning_res, verbose = F) |> as.data.table()
          
          colnames(outl_variog) <- paste0("_", round(distx,0), "_km")
          
          outl_list[["variog"]] <- outl_variog
          
        }
        
      }  
    }
    outliers_results <- do_call(cbind, outl_list)
    colnames(outliers_results) <- paste0("geo.", colnames(outliers_results))
    outliers_results <- ((outliers_results)*1) |> as.data.table()
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