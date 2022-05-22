############ R CODES ######################################################
### Manuscript: Measles epidemic in Southern Vietnam:          ############
### An age-stratified spatio-temporal model for infectious disease counts #
### Authors: Thi Huyen Trang Nguyen     ###################################
### Date: 05 April 2022                 ###################################
###########################################################################


#==================================================================
####################
#### MAIN CODES ####
####################

# LIBRARY
library(surveillance)
library(hhh4contacts)
library(hhh4addon)
library(gridExtra)

library(data.table)
library(ggplot2)

memory.size() ### Checking your memory size
memory.limit() ## Checking the set limit
memory.limit(size=56000) ### expanding your memory

#==================================================================
# DATA IMPORT
measles_sts <- readRDS("measles_sts.rds")       #Import measles data & map & population
contactVN_origin <- readRDS("contactVN.rds")    #source: http://www.socialcontactdata.org/socrates/
pop <- readRDS("populationVN.rds")              #population in VN in 2005, 2009 & 2019

#Population 2005 retrieved from SOCRATES data (http://www.socialcontactdata.org/socrates/)
#POpulation 2009, 2019 retrieved from census data 2009, 2019 
#(http://tongdieutradanso.vn/12-completed-results-of-the-2019-census.html)
#(https://www.gso.gov.vn/du-lieu-va-so-lieu-thong-ke/2019/03/ket-qua-toan-bo-tong-dieu-tra-dan-so-va-nha-o-viet-nam-nam-2009/)

#==================================================================
# DATA & MODEL SET UP

## Modify contact matrix

#Reference: Arregui S, Aleta A, Sanz J, Moreno Y (2018) Projecting social contact matrices 
#to different demographic structures. PLOS Computational Biology 14(12): e1006638. 
#https://doi.org/10.1371/journal.pcbi.1006638
#Method 2: Density correction

contactVN <- matrix(rep(0,16), nrow=4, ncol=4) #Define new contact matrix
for (i in 1:nrow(contactVN_origin)){
  for (j in 1:ncol(contactVN_origin)){
    contactVN[i,j] = contactVN_origin[i,j]*colSums(pop)[4]*pop[j,5]/(pop[j,4]*colSums(pop)[5])
  }
}

colnames(contactVN) <- colnames(contactVN_origin)
rownames(contactVN) <- rownames(contactVN_origin)
#plotC(contactVN)


contactVN_norm <- contactVN / rowSums(contactVN) #the row-normalized version of contact matrix
contactVN_diag <- structure(diag(ncol(contactVN)),dimnames = dimnames(contactVN)) #no-mixing model
contactVN_power_func <- make_powerC(contactVN_norm, truncate = TRUE) # function for power-adjustment of the contact matrix


## Define age-groups and provinces
provinces <- unique(stratum(measles_sts, 1))
nprovinces <- length(provinces)
groups <- unique(stratum(measles_sts, 2))
ngroups <- length(groups)

## Define neighbourhood matrix
mat_neighbourhood <- neighbourhood(measles_sts) #diagonal different than 0 (i.e.=1, for local transmission)

## Define Lunar New Year (2019 & 2020)
newyear <- list(t = epoch(measles_sts) - 1,
                lunar = 1*(epoch(measles_sts) %in% c(398:406,753:759)))
str(newyear)

## Define model matrix by age group indicators
mat_group <- sapply(groups, function (group) {
  list <- which(stratum(measles_sts, which = 2) == group)
  col <- col(measles_sts)
  col[] <- col %in% list
  col}, simplify = FALSE, USE.NAMES = TRUE)
str(mat_group)

## Define model matrix by province indicators
mat_province <- sapply(provinces, function (province) {
  list <- which(stratum(measles_sts, which = 1) == province)
  col <- col(measles_sts)
  col[] <- col %in% list
  col}, simplify = FALSE, USE.NAMES = TRUE)
str(mat_province)


#7 provinces
aa <- mat_province
mat_province1 <- matrix(mapply(sum, aa$`Lam Dong`,aa$`Ho Chi Minh City`,
                               aa$`Binh Phuoc`,aa$`Ba Ria Vung Tau`,
                               aa$`Binh Duong`, aa$`Tay Ninh`, aa$`Dong Nai`,
                               MoreArgs=list(na.rm=T)),ncol=80)  ## provinces in South East area
#13 provinces
mat_province2 <- matrix(mapply(sum, aa$`Long An`, aa$`Dong Thap`, aa$`An Giang`,aa$`Tien Giang`,
                               aa$`Vinh Long`, aa$`Ben Tre`, aa$`Kien Giang`, aa$`Can Tho City`,
                               aa$`Hau Giang`, aa$`Tra Vinh`, aa$`Soc Trang`, aa$`Bac Lieu`,
                               aa$`Ca Mau`,
                               MoreArgs=list(na.rm=T)),ncol=80) ## provinces in Mekong River Delta area
mat_provinceG <- list(SouthEast = mat_province1, Mekong=mat_province2)
str(mat_provinceG)

## Define model matrix by overal sine-cosine effects 
mat_season <- with(c(newyear), unlist(lapply(
  X = 1,
  FUN = function (x) {
    list_indicator <- list(x * sin(2 * pi * t/365),
                           x * cos(2 * pi * t/365))
    names(list_indicator) <- paste0(c("sin", "cos"), "(2*pi*t/365)")
    list_indicator}), recursive = FALSE, use.names = TRUE))
str(mat_season)

## Define groups & province name
nameG <- paste0("`", groups, "`")
nameP <- paste0("`", provinces, "`")
namePG <- paste0("`", names(mat_provinceG), "`")
nameS <- paste0("`", names(mat_season), "`")

## Define Offset
offset_pop <- population(measles_sts)


#==================================================================
#==================================================================
# MODEL FIT
#==================================================================
## Endemic model
f_end <- reformulate(c(1, nameG[-1], namePG[-1], "t","lunar", nameS),
                     intercept = TRUE) 

control<- list(
  end = list(f = f_end,
             offset = offset_pop),
  family = factor(stratum(measles_sts, which = 2)), # group-specific dispersion
  data = c(mat_group, mat_province, newyear, mat_season, mat_provinceG))     
 
fit_endemic <- hhh4(measles_sts, control = control)
summary(fit_endemic)
save(fit_endemic, file="model_fit/fit_endemic.rda")


## Epidemic model (2 components)
f_ne <- reformulate(c(1, nameG[-1], nameP[-1], "log(pop)","lunar"),
                    intercept = TRUE) 

### Fit power law all age groups
fit_PL <- update(fit_endemic,
                 ne = list(f = f_ne,
                           weights = W_powerlaw(maxlag = max(mat_neighbourhood), 
                                                log = TRUE, normalize = FALSE,
                                                initial = c("logd" = log(2))),
                           scale = expandC(contactVN_norm, nprovinces),
                           normalize = TRUE),
                 data = list(pop = population(measles_sts),
                             newyear,mat_group, mat_province)) 

summary(fit_PL)
save(fit_PL, file="model_fit/fit_PL.rda")


### Fit group-specific powerlaw
fit_PL <- get(load("model_fit/fit_PL.rda"))
weights <- addGroups2WFUN(WFUN = fit_PL$control$ne$weights,
                          groups = factor(stratum(measles_sts, which = 2)))
fit_PLG <- update(fit_PL, ne = list(weights = weights))
summary(fit_PLG)
save(fit_PLG, file="model_fit/fit_PLG.rda")

### Fit model with no mixing between groups
fit_nomixing <- update(fit_PL,
                       ne = list(scale = expandC(contactVN_diag, nprovinces)),
                       use.estimates = FALSE)
summary(fit_nomixing)
save(fit_nomixing, file="model_fit/fit_nomixing.rda")

### Fit model with homogeneous mixing
fit_homo <- update(fit_PL,ne = list(scale = NULL),
                   use.estimates = FALSE)
summary(fit_homo)
save(fit_homo, file="model_fit/fit_homo.rda")


### Fit power
fit_PLpower <- fitC(fit_PL, contactVN, 
                    normalize = TRUE, truncate = TRUE)
summary(fit_PLpower)
save(fit_PLpower, file="model_fit/fit_PLpower.rda")

fit_PLGpower <- fitC(fit_PLG, contactVN, 
                    normalize = TRUE, truncate = TRUE)
summary(fit_PLGpower)
save(fit_PLGpower, file="model_fit/fit_PLGpower.rda")


## => final model selection
## Endemic: 1 + group + mekong + trend t + lunar + sincos + offset (pop)
## Epidemic: 1 + group + province + log(pop) + lunar

#==================================================================
##############################
#### SENSITIVITY ANALYSIS ####
##############################


# WEEKLY DATA
#==================================================================
## STS object
measles_stsw <- aggregate(measles_sts, nfreq = 365/7) #aggregate sts object to weekly data
plot(measles_stsw, type = observed ~ time,
     main = "Aggregated over all provinces and age groups")

## Set up
### Define Lunar New Year (2019 & 2020)
newyearw <- list(t = epoch(measles_stsw) - 1,
                 lunar = 1*(epoch(measles_stsw) %in% c(57:58,108:109)))

### Define model matrix by age group indicators
mat_groupw <- sapply(groups, function (group) {
  list <- which(stratum(measles_stsw, which = 2) == group)
  col <- col(measles_stsw)
  col[] <- col %in% list
  col}, simplify = FALSE, USE.NAMES = TRUE)
str(mat_groupw)

### Define model matrix by province indicators
mat_provincew <- sapply(provinces, function (province) {
  list <- which(stratum(measles_stsw, which = 1) == province)
  col <- col(measles_stsw)
  col[] <- col %in% list
  col}, simplify = FALSE, USE.NAMES = TRUE)
str(mat_provincew)

#7 provinces
aa <- mat_provincew
mat_province1 <- matrix(mapply(sum, aa$`Lam Dong`,aa$`Ho Chi Minh City`,
                               aa$`Binh Phuoc`,aa$`Ba Ria Vung Tau`,
                               aa$`Binh Duong`, aa$`Tay Ninh`, aa$`Dong Nai`,
                               MoreArgs=list(na.rm=T)),ncol=80)
#13 provinces
mat_province2 <- matrix(mapply(sum, aa$`Long An`, aa$`Dong Thap`, aa$`An Giang`,aa$`Tien Giang`,
                               aa$`Vinh Long`, aa$`Ben Tre`, aa$`Kien Giang`, aa$`Can Tho City`,
                               aa$`Hau Giang`, aa$`Tra Vinh`, aa$`Soc Trang`, aa$`Bac Lieu`,
                               aa$`Ca Mau`,
                               MoreArgs=list(na.rm=T)),ncol=80)
mat_provinceGw <- list(SouthEast = mat_province1, Mekong=mat_province2)
str(mat_provinceGw)


## Define model matrix by age-specific sine-cosine effects 
mat_seasonw <- with(c(newyearw), unlist(lapply(
  X = 1,
  FUN = function (x) {
    #group_indicator <- get(group)
    list_indicator <- list(x * sin(2 * pi * t/52),
                           x * cos(2 * pi * t/52))
    names(list_indicator) <- paste0(c("sin", "cos"), "(2*pi*t/52)")
    list_indicator}), recursive = FALSE, use.names = TRUE))
str(mat_seasonw)


## Define Offset
offset_popw <- population(measles_stsw) / 100000 #rowSums(population(measles_sts))

## Define groups & province name
nameG <- paste0("`", groups, "`")
nameP = paste0("`", provinces, "`")
namePGw <- paste0("`", names(mat_provinceGw), "`")
nameSw = paste0("`", names(mat_seasonw), "`")
#==================================================================


# DIFFERENT TYPE OF CONTACT MATRICES
#==================================================================

## Original contact matrix 
contactVN_origin
contactVN_origin_norm <- contactVN_origin / rowSums(contactVN_origin) #the row-normalized version of contact matrix

## Contact rate matrix in 2019 & 2009
contact_2019 = contactVN_origin
contact_2019[,1] = contact_2019[,1]/pop[1,3]
contact_2019[,2] = contact_2019[,2]/pop[2,3]
contact_2019[,3] = contact_2019[,3]/pop[3,3]
contact_2019[,4] = contact_2019[,4]/pop[4,3]

contact_2009 = contactVN_origin
contact_2009[,1] = contact_2009[,1]/pop[1,2]
contact_2009[,2] = contact_2009[,2]/pop[2,2]
contact_2009[,3] = contact_2009[,3]/pop[3,2]
contact_2009[,4] = contact_2009[,4]/pop[4,2]

contact_2019norm <- contact_2019 / rowSums(contact_2019) #the row-normalized version of contact matrix
contact_2009norm <- contact_2009 / rowSums(contact_2009) #the row-normalized version of contact matrix



## Results of models fitted in sensitivity analysis can be found in folder "model_fit"
## and retrieved by, for example: fit <- get(load("model_fit/file name of model.rda"))
