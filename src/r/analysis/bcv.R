
# OUTRIDER was switched to branch subsetETraining

devtools::load_all("../OUTRIDER/")
library(ggplot2)
library(LSD)
library(data.table)
library(magrittr)

create_means_cv_dt <- function(ODS, base = FALSE){
    if(isFALSE(base)){
        exp_values <- normalizationFactors(ODS) 
        } else exp_values <- counts(ODS) / normalizationFactors(ODS)
    
    dt <- data.table(
        gene = row.names(ODS),
        theta = theta(ODS),
        geneMeans = rowMeans(exp_values),
        geneMedians = rowMedians(exp_values),
        varMeans = rowMeans(exp_values + (exp_values)^2 / theta(ODS)),
        varMedians = rowMedians(exp_values + (exp_values)^2 / theta(ODS)),
        BCV = 1/sqrt(theta(ODS)))
    return(dt)
}

ods_ss <- readRDS("/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/ss/ods.Rds")
dt_ss <- create_means_cv_dt(ods_ss)
ggplot(dt_ss, aes(geneMeans, varMeans)) + geom_point() + scale_x_log10() + scale_y_log10() + geom_abline()

ggplot(dt_ss, aes(geneMedians, varMedians)) + geom_point() + scale_x_log10() + scale_y_log10() + geom_abline()
ggplot(dt_ss, aes(theta)) + geom_histogram()
ggplot(dt_ss, aes(1/sqrt(theta))) + geom_histogram()

ggplot(dt_ss, aes(BCV, geneMeans)) + geom_point() + scale_y_log10()

ggplot(dt_ss, aes(geneMeans, BCV)) + geom_point() + scale_x_log10() + ylim(c(0,0.5)) + theme_bw()
ggplot(dt_ss, aes(geneMedians, BCV)) + geom_hex() + scale_x_log10() + ylim(c(0,0.5)) + theme_bw() + ggtitle("SS all")


# Do the same with counts normalized using size factors only
ods_base <- readRDS("/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/ss/ods_unfitted.Rds")
ods_base <- ods_base[mcols(ods_base)$passedFilter, ]
ods_base <- estimateSizeFactors(ods_base)
row.names(ods_base) <- rowData(ods_base)$gene_name_unique
normalizationFactors(ods_base) <- matrix(sizeFactors(ods_base), nrow = nrow(ods_base), ncol=ncol(ods_base), byrow = TRUE)
# ods_base <- fit(ods_base)   # Do in ouga
saveRDS(ods_base, "/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/ss/ods_base.Rds")
ods_base <- readRDS("/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/ss/ods_base.Rds")
dt_base <- create_means_cv_dt(ods_base, base = FALSE)

ggplot(dt_base, aes(geneMedians, BCV)) + geom_hex() + scale_x_log10() # + ylim(c(0,1))
heatscatter(dt_ss$BCV, dt_base$BCV); grid(); abline(0,1)
boxplot(theta(ods_ss), theta(ods_base), outline = F)



###### Artery Aorta #####
ods_artery_aorta <- readRDS('/s/project/scared/paper/revision/run0910-final/data/fitOutrider/NLas_TCNo/Artery_Aorta_ODS.RDS')
dim(ods_artery_aorta)   # 224 samples
dt_artery_aorta <- create_means_cv_dt(ods_artery_aorta)
ggplot(dt_artery_aorta, aes(geneMedians, BCV)) + geom_hex() + scale_x_log10() + ylim(c(0,0.5))
normalizationFactors(ods_artery_aorta) <- t(t(counts(ods_artery_aorta)+1)/sizeFactors(ods_artery_aorta))
# ods_artery_aorta <- fit(ods_artery_aorta)   # Do in ouga
ods_artery_aorta_base<- readRDS('/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/ss/ods_artery_aorta.Rds')
dt_artery_aorta_base <- create_means_cv_dt(ods_artery_aorta_base)
ggplot(dt_artery_aorta_base, aes(geneMedians, BCV)) + geom_hex() + scale_x_log10() + ylim(c(0,0.5))


###### Brain anterior #####
ods_brain_anterior <- readRDS('/s/project/scared/paper/revision/run0910-final/data/fitOutrider/NLas_TCNo/Brain_Anterior_cingulate_cortex_BA24_ODS.RDS')
dim(ods_brain_anterior)   # 84 samples
dt_brain_anterior <- create_means_cv_dt(ods_brain_anterior)
# ggplot(dt_brain_anterior, aes(geneMedians, BCV)) + geom_hex() + scale_x_log10() + ylim(c(0,0.5))
normalizationFactors(ods_brain_anterior) <- t(t(counts(ods_brain_anterior)+1)/sizeFactors(ods_brain_anterior))
# ods_brain_anterior <- fit(ods_brain_anterior)   # Do in ouga
ods_brain_anterior_base<- readRDS('/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/ss/ods_brain_anterior.Rds')
dt_brain_anterior_base <- create_means_cv_dt(ods_brain_anterior_base)
ggplot(dt_brain_anterior_base, aes(geneMedians, BCV)) + geom_hex() + scale_x_log10() # + ylim(c(0,0.5))




ods <- makeExampleOutriderDataSet(n = 2000)
dim(ods)
ods_test <- ods_ss

## Execute in ouga
ods_list <- list()
ods_list <- bplapply(1:5, function(i){
    set.seed(i)
    ods_test <- ods_ss
    genes_include <- sample(c(T,F), replace = T, prob = c(.1,.9), nrow(ods_test))
    featureExclusionMask(ods_test) <- genes_include
    ods_test <- OUTRIDER(ods_test)
    ods_test
    }, 
    BPPARAM = MulticoreParam(workers = 20)
)

ods_list <- readRDS('/s/project/genetic_diagnosis/tmp_yepez/ods_list_5simulations.Rds')

dt_l <- lapply(ods_list, create_means_cv_dt)
names(dt_l) <- paste0("sim_", 1:5)
dt_sim_full <- rbindlist(dt_l, idcol = "sim")
dt_ss[, sim := 'or']
dt_sim_full <- rbind(dt_ss, dt_sim_full)

ggplot(dt_sim_full, aes(geneMedians, BCV)) + geom_hex() + scale_x_log10() + ylim(c(0,0.5)) + facet_wrap(~sim) + scale_fill_gradient(trans = 'log10')

res_l <- lapply(ods_list, OUTRIDER::results)
res_sim_full <- rbindlist(res_l, idcol = "sim")
res_ss <- OUTRIDER::results(ods_ss)
res_ss[, sim := 'or']
res_sim_full <- rbind(res_ss, res_sim_full)
res_sim_full[, .N, by = sim]
setdiff(res_sim_full[padjust < .001 & sim == '5', geneID], res_sim_full[sim == 'or', geneID])
setdiff(res_sim_full[padjust < .001 & sim == 'or', geneID], res_sim_full[sim == '1', geneID])


dim(E(ods_test))
dim(D(ods_test))
E(ods_test)

dc <- dcast(dt_sim_full[,.(gene, theta, sim)], formula = gene ~ sim, value.var = 'theta')

heatpairs(as.matrix(dc[, 2:7]), log = 'xy') 
