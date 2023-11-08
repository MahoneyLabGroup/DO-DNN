library(qtl2)
library(parallel)
detectCores()
setwd("/projects/korstanje-lab/chenm/DO-DNN/Churchill-046_JAC_DO_Aging-MegaMUGA")

cross <- read_cross2("Churchill-046_JAC_DO_Aging-MegaMUGA.json")
pr <- calc_genoprob(cross , error_prob=0.002, cores = detectCores())
#saveRDS(pr, file = "pr.rds")
#pr<-readRDS(file = "pr.rds")


apr <- genoprob_to_alleleprob(pr, cores = detectCores())

k <- calc_kinship(apr, "loco")
covar = read.csv("/projects/korstanje-lab/chenm/DO-DNN/Churchill-046_JAC_DO_Aging-MegaMUGA/Churchill-046_JAC_DO_Aging-MegaMUGA_covarbatchMME.csv")

get_dat_dummpy<-function(dat, name){
  dat_unique = unique(dat)
  num_column = dim(dat_unique)[1] - 1
  num_row = dim(dat)[1]
  dat_dummy = matrix(0, nrow = num_row, ncol = num_column)
  print(dat_unique)
  col_names = c()
  for (i in 1:num_column) {
    d = dat_unique[[1]][i]
    print(d)
    coln = paste0(name, d)
    dat_dummy[dat == d,i] = 1
    col_names = append(col_names, coln)
  }
  colnames(dat_dummy) = col_names
  rownames(dat_dummy) = rownames(dat)
  return(dat_dummy)
}


df_ave_3 = read.csv( "/projects/korstanje-lab/chenm/DO-DNN/data_zarrs_3/df_ave_3_log1698952711.csv")
meta <- merge(covar[c("AnimalID", "Sex", "gen", "Location", "AgeNum")],df_ave_3,by.x = "AnimalID", by.y="Animal_ID")

###NOTE: DO-1224  HAS TWO IMAGES: 334414_DO_1224_B_PAS (2).svs [0].zarr AND 334235_DO_1224_B_PAS.svs [0].zarr
meta[meta$AnimalID == "DO-1224",]
meta = meta[!meta$filename == "334235_DO_1224_B_PAS.svs [0].zarr",]

### generate covariance matrix
sex <- meta[, "Sex", drop=FALSE]
rownames(sex) =  meta[,"AnimalID"]

gen = meta[,"gen", drop=FALSE]
rownames(gen) = meta[,"AnimalID"]

age = age = meta[,"age", drop=FALSE]
rownames(age) = meta[,"AnimalID"]

location = meta[,"Location", drop=FALSE]
rownames(location) = meta[,"AnimalID"]

dataset = meta[,"dataset", drop=FALSE]
rownames(dataset) = meta[,"AnimalID"]


sex = get_dat_dummpy(sex, "sex")
gen = get_dat_dummpy(gen, "gen")
age = get_dat_dummpy(age, "age")
location = get_dat_dummpy(location, "location")
dataset = get_dat_dummpy(dataset, "dataset")

covar = sex
covar = cbind(covar, gen)
covar = cbind(covar, age)
covar = cbind(covar, location)
covar = cbind(covar, dataset)


score = meta[,"pc1", drop=FALSE]
rownames(score) = meta[,"AnimalID"]

out <- scan1(apr, score, k, addcovar = covar, cores= detectCores())

operm <- scan1perm(apr, score,k,  addcovar = covar, n_perm=100, cores=detectCores())

