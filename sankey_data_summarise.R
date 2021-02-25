library(data.table)
library(dplyr)

edgeDat = fread("shankey_edge_data.csv") %>% as.data.frame()
edgeDat.2 = edgeDat %>% dplyr::filter(group == "dr_pathway" & 
                                        target != "TGF-BETA SIGNALING")
edgeDat.3 = dplyr::inner_join(edgeDat.2, edgeDat, by = c("target" = "source"))
# edgeDat.4 = edgeDat.3 %>% 
#   dplyr::group_by(source, target) %>% 
#   dplyr::summarise_at(.vars = vars("target.y"), .funs = paste0(collapse = ", "))

edgeDat.4 = stats::aggregate(target.y ~., edgeDat.3, toString)
edgeDat.5 = stats::aggregate(source ~., edgeDat.4, toString) %>% 
  dplyr::select(c(5,1,4))
colnames(edgeDat.5) = c("TGF-β driver genes", "Other perturbed STPs", "Enriched Breast Cancer genes")
fwrite(edgeDat.5, file="global_perturbation_summarised.csv")
