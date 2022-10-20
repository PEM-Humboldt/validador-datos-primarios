# B. Modulo ambiental

# B.2 etiquetado ambiental
do_environmental_tag <- function(
    env = envars, data_base, col_lon, col_lat, univar = F, multivar = F,
    test_univar = c("zscore", "std", "iqr", "rjack", "mad" ), test_multivar = c("pca_error", "maha"),
    univar_details = list("mtpl_std" = 3, "mtpl_iqr" = c(1.5, 3, 5), "mtpl_mad" = c(1.5, 3, 5)), 
    multivar_details = c("pval" = 0.95), min_occs = 7
){
  
  env_space <- gen_env_space(env. = env, data.base = data_base, col.lon = col_lon, col.lat = col_lat)
  
  if(univar == T){
    
    # inicializar output
    
    outl_list <- list()
    
    if(nrow(data_base) <=  min_occs){
      return(out_univar)
      
    }else{
      
      for(i in 1:ncol(env_space)){
        
        x <-  env_space[, ..i] |> as.matrix() |> as.vector()
        xi <- na.omit(x)
        xNA <- which(is.na(x) == T)
        
        out_ <- list()
        
        for(a in 1:length(test_univar)){
          
          test_a <- test_univar[a]
          
          if (test_a == "zscore"){
            z <- calc_zscore(x = xi )
            if(length(xNA)!= 0){
              for(f in 1:length(xNA)){
                z <- append(x = z, values = NA, after = xNA[f]-1 )  
              }
            }
            out_[["zscore"]] <- z
          } 
          
          if (test_a == "std"){
            mltpl <- univar_details[["mtpl_std"]]
            std_v <- ((!calc_std(x = xi , threshold = mltpl))*1)
            if(length(xNA)!= 0){
              for(f in 1:length(xNA)){
                std_v <- append(x = std_v, values = NA, after = xNA[f]-1 )  
              }
            }
            out_[["std"]] <- std_v
          }
          
          if (test_a == "rjack"){
            rjackv <- (!calc_rjack(x = xi ))*1
            if(length(xNA)!= 0){
              for(f in 1:length(xNA)){
                rjackv <- append(x = rjackv, values = NA, after = xNA[f]-1 ) 
              }
            }
            out_[["rjack"]] <- rjackv
          }
          
          if(test_a == "iqr"){
            mltpl <- univar_details[["mtpl_iqr"]]
            out_a <- list()
            for(c in 1:length(mltpl)){
              out_c <- (!iqr(x = xi, mtp = mltpl[c]))*1
              if(length(xNA)!= 0){
                for(f in 1:length(xNA)){
                  out_c <- append(x = out_c, values = NA, after = xNA[f]-1 )
                }
              }
              out_a[[c]] <- out_c
            }
            out_a <- do_call(cbind, out_a) |> as.data.table()
            colnames(out_a) <- paste0("_", mltpl)
            
            out_[["iqr"]] <- out_a
          }
          
          if(test_a == "mad"){
            mltpl <- univar_details[["mtpl_mad"]]
            out_a <- list()
            for(c in 1:length(mltpl)){
              out_c <- (!calc_mad(x = xi, threshold = mltpl[c]))*1 
              if(length(xNA)!= 0){
                for(f in 1:length(xNA)){
                  out_c <- append(x = out_c, values = NA, after = xNA[f]-1 )
                }
              }
              out_a[[c]] <- out_c
            }
            out_a <- do_call(cbind, out_a) |> as.data.table()
            colnames(out_a) <- paste0("_", mltpl)
            
            out_[["mad"]] <- out_a
          }
        }
        out_ <- do_call(cbind, out_)
        var <- colnames(env_space)[i]
        outl_list[[var]] <- out_
      }
    }
    
    univar_results <- do_call(cbind, outl_list)
  }
  
  if(multivar == T){
    
    out_multivar <- list()
    
    # hallar los datos perdidos, no pueden ser trabajados en los test multivariados
    xNA <- apply(X = env_space, 2, FUN = function(X){which(is.na(X))}) |> unlist() |> unique()
    
    if(is.matrix(xNA)){
      xNA <- xNA %>% as.vector() %>% unique()
    }
    
    # media y escalar las variables
    if(length(xNA)!= 0){
      mu <-  colMeans(env_space[-xNA, ])
      env_space_scale <-  scale(env_space[-xNA, ], center = mu, scale = T)  
    }else{
      mu <-  colMeans(env_space)
      env_space_scale <-  scale(env_space, center = mu, scale = T)
    }
    
    
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
      if(length(xNA)!= 0){
        for(f in 1:length(xNA)){
          errorq_pca <- append(x = errorq_pca, values = NA, after = xNA[f]-1 )
        }
      }
      
      out_multivar[["pca_error"]] <- errorq_pca
      
    }
    
    if("maha" %in% test_multivar){
      
      if(!require(ClassDiscovery)) install.packages("ClassDiscovery")
      
      
      # establecer la distancia de mahalanobis y hallar valores de probabilidad
      dist.maha <- mahalanobis(env_space_scale, colMeans(env_space_scale), cov(env_space_scale)) %>% round(2)
      pvalue <- pchisq(dist.maha, df=ncol(env_space_scale)-1, lower.tail=FALSE)
      
      # binarizar los valores de probabilidad, marcando los que esten por debajo
      # de significatividad estadistica deseada 1 - pvalor
      
      pvalue <- !(pvalue <= (1- multivar_details["pval"]))
      pvalue <- pvalue*1
      
      # agregar valores NA
      if(length(xNA)!= 0){
        for(f in 1:length(xNA)){
          dist.maha <- append(x = dist.maha, values = NA, after = xNA[f]-1 )
          pvalue <- append(x = pvalue, values = NA, after = xNA[f]-1 )
        }
      }
      
      maha <- cbind(dist.maha, pvalue) |> as.data.frame()
      
      # adjuntar a la data.frame de salida
      out_multivar[["maha"]] <- maha$dist.maha
      out_multivar[["pvalue"]] <- maha$pvalue
    }
    
    multivar_results <- do_call(cbind, out_multivar)
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