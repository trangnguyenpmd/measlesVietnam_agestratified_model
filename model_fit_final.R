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
#plotC(contactVN)  #plot the contact matrix


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
for (i in 1:4){
  contact_2019[,i] = contact_2019[,i]/pop[i,3]
}

contact_2009 = contactVN_origin
for (i in 1:4){
  contact_2009[,i] = contact_2009[,i]/pop[i,2]
}

contact_2019norm <- contact_2019 / rowSums(contact_2019) #the row-normalized version of contact matrix
contact_2009norm <- contact_2009 / rowSums(contact_2009) #the row-normalized version of contact matrix



## Results of models fitted in sensitivity analysis can be found in folder "model_fit"
## and retrieved by, for example: fit <- get(load("model_fit/file name of model.rda"))


#==================================================================
##############################
#### EXTRACT COEFFICIENTS #############
##############################

fit <- get(load("model_fit/fit_PLG.rda"))
fit$dim
summary(fit)
summary(fit, idx2Exp = TRUE)

nterms <- terms(fit)$nGroups #+ 8
coefs <- coef(fit)[1:nterms]
CIs <- confint(fit)[1:nterms, ]
id_log <-  c(grep("over", names(coefs)), grep("neweights.d", names(coefs)))
coefs[id_log] <- log(coefs[id_log])
CIs[id_log, ] <- log(CIs[id_log, ])
tab <- round(cbind(coefs, CIs), 3)
tab
# write.csv(as.data.frame(tab),"fit.csv", row.names = TRUE)


#==================================================================
##############################
#### MAIN FIGURE #############
##############################

# FUNCTION
# Reference: Bracher, J. and L. Held, Endemic-epidemic models with 
# discrete-time serial interval distributions for infectious disease prediction. 
# International Journal of Forecasting, 2020
# DOI: https://doi.org/10.1016/j.ijforecast.2020.07.002

decompose_epidemic_component <- function(fit){
  # extract info:
  sts <- fit$stsObj
  max_lag <- if(class(fit)[1] == "hhh4lag") fit$max_lag else 1
  subset <- fit$control$subset
  n_units <- ncol(sts@observed)
  param <- hhh4addon:::lambda_tilde_complex_neighbourhood(fit, periodic = FALSE,
                                                          subset = 1:max(subset))
  
  # initialize:
  contributions <- array(dim = c(max(subset),
                                 n_units,
                                 n_units,
                                 max_lag),
                         dimnames = list(1:max(subset),
                                         paste0("from.", colnames(sts@observed)),
                                         paste0("to.", colnames(sts@observed)),
                                         paste0("lag", 1:max_lag)))
  # fill:
  for(t in subset){
    phi_temp <- param$lambda[,,t]
    obs_preceding <- sts@observed[t - max_lag:1, , drop = FALSE]
    for(lag in 1:max_lag){
      inds <- seq(to = n_units*(max_lag - lag + 1), length.out = n_units)
      phi_this_lag <- phi_temp[, inds]
      contributions[t, , , lag] <- t(phi_this_lag)*matrix(obs_preceding[max_lag - lag + 1, ], ncol = n_units, nrow = n_units)
    }
  }
  return(contributions)
}




#----- Write function to construct df based on age groups
#   - 1: 00-04
#   - 2: 05-14
#   - 3: 15-24
#   - 4: 25+
library(tidyverse)

constructDf <- function(agegr) {
  beg            <- 20*(agegr - 1) + 1
  end            <- 20*agegr
  columns        <- c(beg:end)
  fitg           <- as.data.frame(fit_component[[agegr]])
  fitg$observed  <- rowSums(as.data.frame(fit$stsObj@observed[2:912, columns]))  #for weekly, change 912 to 131
  fitg$observed  <- ifelse(fitg$observed == 0, NA, fitg$observed)
  age_gr         <- c("0004", "0514", "1524", "2500")
  age_index      <- age_gr[agegr]
  age_id         <- paste0("from", age_index, ".", "to", age_index)
  fitg$within    <- fit_epi[[age_id]]
  # Create variable "from*" for different age groups
  for (i in 1:4) {
    if (i != agegr) {
      a <- paste0("from", age_gr[i])
      b <- paste0("from", age_gr[i], ".to", age_gr[agegr])
      eval(call("<-", as.name(a), fit_epi[[b]]))
    } else {
      a <- paste0("from", age_gr[agegr])
      eval(call("<-", as.name(a), rep(0, 130)))
    }
  }
  fitg$from0004  <- from0004
  fitg$from0514  <- from0514
  fitg$from1524  <- from1524
  fitg$from2500  <- from2500
  fitg$withina   <- fitg$endemic + fitg$within
  fitg$from0004a <- fitg$withina + fitg$from0004
  fitg$from0514a <- fitg$from0004a + fitg$from0514
  fitg$from1524a <- fitg$from0514a + fitg$from1524
  fitg$from2500a <- fitg$from1524a + fitg$from2500
  return(fitg)
}

#==================================================================
#==================================================================
# Figure 4 (main text): Fitted each of 4 age groups of all provinces
#==================================================================

## Set up data frame for plot in ggplot2

fit <- get(load("model_fit/fit_PLG.rda")) #select model fit from folder "model_fit"

## Extract all component
plotHHH4_fitted_groups(fit,
                       groups = stratum(measles_sts, which = 2), 
                       col = c("grey85","#0072B2","#E69F00"),
                       units = NULL, pch = 20, legend = 2,
                       start = fit$stsObj@start,
                       legend.args = list(legend = c("From other age groups",
                                                     "Within age group", 
                                                     "Endemic"))) -> fit_component

## Extract WHO ACQUIRED INFECTION FROM WHOM
fit_epi <- decompose_epidemic_component(fit)
fit_epi <- fit_epi[1:912,1:80,1:80,]
#fit_epi <- fit_epi[1:131,1:80,1:80,]  # for model fit with weekly data
fit_epi[is.na(fit_epi)] = 0
str(fit_epi)

a1 <- c(rep("from0004",20),rep("from0514",20),rep("from1524",20), rep("from2500",20))
a2 <- c(rep("to0004",20),rep("to0514",20),rep("to1524",20), rep("to2500",20))

dimnames(fit_epi) <- list(c(1:912),a1,a2)
#dimnames(fit_epi) <- list(c(1:131),a1,a2)  # for model fit with weekly data
fit_epi <- setDT(data.frame(fit_epi))
colnames(fit_epi) <- substr(colnames(fit_epi),1,15)
fit_epi <- t(apply(fit_epi,1, function(x) tapply(x,colnames(fit_epi),sum)))
fit_epi <- as.data.frame(fit_epi[2:912,]) 
#fit_epi <- as.data.frame(fit_epi[2:131,])    # for model fit with weekly data


## Create data frame
#----- Apply function above to create 4 dfs
fitg1 <- constructDf(agegr = 1)
fitg2 <- constructDf(agegr = 2)
fitg3 <- constructDf(agegr = 3)
fitg4 <- constructDf(agegr = 4)


#-------------------------------------------------------------------------------
#----- Prepare data for plot
date912 = seq(as.Date("2018-01-01"), as.Date("2020-06-30"), by = "day")
#dateweek = seq(as.Date("2018-01-01"), as.Date("2020-06-30"), by = 7) #for model with weekly data

fitg1$agegr <- "00-04"
fitg2$agegr <- "05-14"
fitg3$agegr <- "15-24"
fitg4$agegr <- "25+"

fitg_all <- rbind(fitg1, fitg2, fitg3, fitg4)
fitg_all$date912 <- rep(date912[2:912], 4)
#fitg_all$dateweek <- rep(dateweek, 4)   #for model with weekly data

head(fitg_all)

#----- Plot
fitg_all %>% select(date912, agegr, from2500a, from1524a, from0514a,
                    from0004a, withina, endemic, observed) %>%
  gather(-c(date912, agegr, observed), key = "indicator", value = "value") %>%
  mutate(indicator = factor(indicator, levels = c("from2500a", "from1524a", "from0514a",
                                                  "from0004a", "withina", "endemic"))) %>%
  ggplot() +
  geom_area(aes(x = date912, y = value, fill = indicator), 
            position = position_dodge(width = 0.001)) +
  geom_point(aes(x = date912, y = observed), color = "black", fill = "black", size = 0.5, shape = 20) +
  facet_wrap(~agegr, scales = "free", ncol = 2) +
  scale_fill_manual(name = "",
                    values = c("from2500a" = "#FFFF00",
                               "from1524a" = "#E69F00", 
                               "from0514a" = "#FF6633", 
                               "from0004a" = "#CC3300", 
                               "withina"   = "#0072B2",
                               "endemic"   = "grey85",
                               "Observed"  = "black"),
                    labels = c("From age group 25+", 
                               "From age group 15-24", 
                               "From age group 05-14", 
                               "From age group 00-04", 
                               "Within age group",
                               "Endemic",
                               "Observed")) +
  scale_x_date(limits = as.Date(c("2017-12-31", "2020-07-31")),
               date_breaks = "1 year",
               date_minor_breaks = "1 month",
               date_labels = "%Y") +
  labs(
    x = "Date of onset",
    y = "No. cases"
  ) +
  theme_minimal() +
  theme(
    axis.line         = element_line(size=0.2, color = "black"),
    axis.text         = element_text(color = "black", size = 10),
    axis.title        = element_text(color = "black", size = 12),
    axis.ticks        = element_line(size = 0.1),
    axis.ticks.length = unit(.05, "cm"),
    legend.position   = c(0.9, 0.9),
    legend.text       = element_text(size = 10), 
    legend.key.size   = unit(0.15, "in"),
    strip.text        = element_text(size = 12, face = "bold")
  ) ->> mainFigure

#----- Save figure
# Note: with tiff format, we should use "compression" to reduce storage space
# otherwise, with publication quality, you will not be able to upload figure to their systems
tiff("mainFigure.tiff", units = "in", width = 10, height = 8, res = 600, compression = "lzw")
mainFigure
dev.off()


#==================================================================
#==================================================================
# Supplementary Figure 1: Fitted all for 4 age groups
#==================================================================

fit <- get(load("model_fit/fit_PLG.rda"))

tiff("Supplementary Figure S1.tiff", units = "in",
     width = 9, height = 7.5, res = 300)

plotHHH4_fitted_groups(fit,
                       groups = stratum(measles_sts, which = 2), total = TRUE, pch = 20,
                       col = c("grey85","#0072B2","#E69F00"),
                       legend.args = list(legend = c("From other age groups", "Within age group", "Endemic"))) -> fit_all

dev.off()


#==================================================================
#==================================================================
# Supplementary Figure 2: Fitted by provinces
#==================================================================

fit <- get(load("model_fit/fit_PLG.rda"))

tiff("Supplementary Figure S2.tiff", units = "in",
     width = 10, height = 8, res = 600)

plotHHH4_fitted_groups(fit,
                       groups = factor(stratum(measles_sts, which = 1), levels = provinces),
                       col = c("grey85","#0072B2","#E69F00"),
                       names = provinces, units = NULL, 
                       legend = 1, 
                       legend.args = 
                         list(cex = 1,
                              legend = c("From other provinces", "Within provinces", "Endemic")),
                       pch = 20) #-> aa

dev.off()
