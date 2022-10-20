# Correlaciones ambientales
do_corr_envars <- function(envars = envars, sample_size = 10000){
  if(!require(terra)) install.packages("terra")
  sample_envars <- terra::spatSample(envars, size = 10000, method = "random", replace = F, na.rm = T,
                                     as.df= T)
  correlacion <- sample_envars |> cor()
  return(correlacion)
}