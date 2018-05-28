#' Gene set for neuronal marker
#' @description  This is a sample example for neuronal gene marker in different category.
#' @docType data
#' @usage data(marker_gene)
#' @format A datset of neuronal gene listed in different category (List).
"marker_gene"

# Data generation:
# genr              <- c("Slc6a1", "Lhx6", "Gad1")
# genr_neuro        <- c("Cdh13","Slc24a2","Gda","Lphn2","Bcl11b","Calb1","Atp1b1","Rgs4", "Pde1a","Penk")
# genr_neuro_pyra1  <- c("Nrn1", "Pcp4", "Rprm", "Enpp2", "Rorb","Rasgrf2", "Wfs1","Fos","Plcxd2","Crym", "X3110035E14Rik","Pvrl3")
# genr_neuro_pyra2  <- c("Cux2", "Kcnk2", "Nr4a2")
# genr_neuro_inter1 <- c("Sst", "Chodl", "Npy", "Reln", "Satb1","Grin3a", "Cort", "Crhbp", "Th","Synpr")
# genr_neuro_inter2 <- c("Pvalb", "Chrm2", "Tac1","Thsd7a")
# genr_neuro_inter3 <- c("Ndnf", "Hapln1", "Kit","Gabrd","Rab3c","Nos1","Lamp5")
# genr_neuro_inter4 <- c("Cxcl14", "Cpne5", "Rgs12")
# genr_neuro_inter5 <- c("Cck","Cnr1","Sema3c","Fxyd6","Slc17a8","Sncg","Rgs10","Nov")
# genr_neuro_inter6 <- c("Htr3a","Vip", "Calb2","Npy2r", "Crh","Pthlh", "Tac2","Qrfpr")
# genr_nonneuro     <- c("Pax6","Plp1","Aldoc","Sulf2", "Serpini1","Zcchc12", "Snca","Id2","Scg2","Gap43","Neurod6","Fam19a1","Enc1", "Arpp21","Calm2","Cox6a2")
# marker_gene <- list(genr = genr, genr_neuro = genr_neuro, genr_neuro_pyra1 = genr_neuro_pyra1, genr_neuro_pyra2 = genr_neuro_pyra2, genr_neuro_inter1 = genr_neuro_inter1, genr_neuro_inter2 = genr_neuro_inter2, genr_neuro_inter3 = genr_neuro_inter3, genr_neuro_inter4 = genr_neuro_inter4, genr_neuro_inter5 = genr_neuro_inter5, genr_neuro_inter6 = genr_neuro_inter6, genr_nonneuro = genr_nonneuro)
# devtools::use_data(marker_gene,overwrite = TRUE,compress = "xz")


#' Single cell RNA seq data
#' @description  Single cell RNA seq data for mouse hippocampus
#' @docType data
#' @usage data(single_cell)
#' @format A datset of single cell RNA seq data.
"single_cell"
# devtools::use_data(single_cell,overwrite = TRUE,compress = "xz")
