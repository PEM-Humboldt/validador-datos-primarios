# Este script busca extraer del conjunto de datos denominado set16 (version 2021) aquellas
# especies con registros curados por expertos con el fin de desarrollar una aplicación demo
# del protocolo de limpieza de datos.

library(data.table)
library(sf)

# datos curados por expertos
registros_expertos_noXenartos <- data.table::fread("registros_demo/registros_expertos_noXenartos.csv", 
                                       stringsAsFactors = F)
spp_expertos <- unique(registros_expertos_noXenartos$acceptedNameUsage)
# 1948 spp
# datos sin curar (set16)

# load("registros_demo/set16.RData")
# set16 <- set16 |> as.data.table()
# subset16 <- set16[acceptedNameUsage %in% spp_expertos]

#rm(set16); gc()

spp_subset16 <- unique(subset16$acceptedNameUsage)
#1871 spp

# diferencia 77 especies

spp_expertos[!spp_expertos %in% spp_subset16]

# [1] "Trachemys medemi"                 "Cavendishia dendrophila"         
# [3] "Thibaudia aurantia"               "Cavendishia sophoclesioides"     
# [5] "Orthaea caudata"                  "Gomphichis cladotricha"          
# [7] "Lepanthes schnitteri"             "Liparis schneideri"              
# [9] "Malaxis mucronulata"              "Masdevallia aenigma"             
# [11] "Masdevallia sumapazensis"         "Microchilus quetamensis"         
# [13] "Restrepia pandurata"              "Stelis cycloglossa"              
# [15] "Stelis umbriae"                   "Telipogon alvarezii"             
# [17] "Telipogon roseus"                 "Trichosalpinx escobarii"         
# [19] "Diaemus youngii"                  "Molossops temmincki"             
# [21] "Neoplatymops mattogrossensis"     "Vampyrodes major"                
# [23] "Balantiopteryx infusca"           "Myotis caucensis"                
# [25] "Lichonycteris degener"            "Coccothrinax argentata"          
# [27] "Landoltia punctata"               "Acoelorrhaphe wrightii"          
# [29] "Saucerottia castaneiventris"      "Orthemis discolor"               
# [31] "Ischnura chingaza"                "Ischnura cruzi"                  
# [33] "Mesamphiagrion risi"              "Atelopus marinkellei"            
# [35] "Magnolia caricaefragrans"         "Magnolia cararensis"             
# [37] "Magnolia neomagnifolia"           "Boana picturata"                 
# [39] "Pitcairnia huilensis"             "Berberis meollancensis"          
# [41] "Espeletia pisbana"                "Fernandezia ortiziana"           
# [43] "Espeletia mutabilis"              "Espeletiopsis diazii"            
# [45] "Espeletia tapirophila"            "Espeletiopsis trianae"           
# [47] "Maxillaria pastoense"             "Crossoglossa exigua"             
# [49] "Aspidogyne carauchana"            "Masdevallia renzii"              
# [51] "Pleurothallis foliosa"            "Stelis pastoensis"               
# [53] "Ophidion carrilloi"               "Temporary species"               
# [55] "Pleurothallis sibatensis"         "Bolitoglossa tamaense"           
# [57] "Atelopus tamaense"                "Blakea vallensis"                
# [59] "Licaria colombiana"               "Meriania trianae"                
# [61] "Ocotea celastroides"              "Ocotea cuatrecasasii"            
# [63] "Oreopanax niger"                  "Oreopanax velutinus"             
# [65] "Schefflera bifurcata"             "Swartzia amabale"                
# [67] "Swartzia radiale"                 "Williamodendron quadrilocellatum"
# [69] "Bdelyrus laplanadae"              "Cryptocanthon buriticaorum"      
# [71] "Coendou ichillus"                 "Coendou vestitus"                
# [73] "Myoprocta acouchi"                "Eriocnemis isabellae"            
# [75] "Cerdocyosi thous"                 "Cebus leucocephalus"             
# [77] "Aotus grisemembra" 

save(subset16, file = "registros_demo/subset16.RData", row.names = F)
