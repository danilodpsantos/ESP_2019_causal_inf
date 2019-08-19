setwd('U:/ESP2019/ESP48/AdditionalMaterial_R/dat') #insert folder path where data are saved locally
dat <- read.csv('lab3.csv')

# Get an idea of the data set
summary(dat)
nrow(dat)
nrow(dat[!duplicated(dat[,'id']),])
subset(dat,id==93003)[,c('id','month','art','cd4_v','min_cd4_v')]
subset(dat,id==93006)[,c('id','month','art','cd4_v','min_cd4_v')]

###############################
# Here: Exercise Question 1-3 #
###############################

# 1) 
length(unique(dat$id))
nrow(dat[!duplicated(dat[,'id']),])

# 2) Patient 93003 starts ART in month 67 when CD4 cell count is 732. 
# Note also at month 0 the patient's CD4 cell count is already below 500 and in month 49 it is below 450.
# Therefore, the patient never follows the treatment strategy for the 500 threshold, follows the treatment strategy
# for 450 for 49 months (but fails to start in month 49 as appropriate for this strategy), 
# and follows the other treatment strategies for 67 months (but starts prematurely in month 67 for these strategies).
subset(dat,id==93003)[,c('id','month','art','cd4_v','min_cd4_v')]

# 3) Patient 93006 starts ART right away in month 0 when CD4 count is 310. Therefore, this patient is following 
# strategies for the 500, 450, 400, and 350 threshold, but started too early to follow 300, 250, or 200. 
subset(dat,id==93006)[,c('id','month','art','cd4_v','min_cd4_v')]

#############################
# Let's continue analysis   #
#############################

# IPTW denominator model (on the subset of those months were eligible for ART, i.e. month >=0 ???)
mod <- glm(art ~ month + monthsq
	+ factor(age_0_cat) + SEX + factor(origin) + factor(mode)
	+ factor(year_0_cat) + factor(cd4_0_cat) + factor(rna_0_cat)
	+ factor(cd4_v_cat) + factor(rna_v_cat),
	family = binomial(), data = subset(dat,eligible2==1))

dat$probA.d <- ifelse(dat$eligible2==1,	predict(mod, type = 'response'),NA)

# Question 4
summary(dat$probA.d)

# Clone strategy
for(i in 1:7) {
	dat[[paste('regime',i,sep='')]] <- i
}

head(dat)

datrep <- reshape(dat,varying=c('regime1','regime2','regime3','regime4','regime5','regime6','regime7'),
                  v.names='regime', idvar=c('id','month'),	direction='long')

head(datrep)

# helper variables
datrep$threshold <- 200+(datrep$regime-1)*50
datrep$newid <- with(datrep,interaction(id,regime))

summary(factor(subset(datrep,month==0)$threshold))
summary(factor(datrep$threshold))

# IMPORTANT part: remove clones who are not following a respective treatment strategy at BASELINE
# code is a bit complex here 
datrep$min.allowed <- ifelse(datrep$art_base==1,50*ceiling((datrep$cd4_0)/50),0)
datrep <- subset(datrep, threshold>=min.allowed)
datrep$max.allowed <- ifelse(datrep$art_base==0, 50*floor(datrep$cd4_0/50),500)
datrep <- subset(datrep, threshold<=max.allowed)

summary(factor(subset(datrep,month==0)$threshold))

# Example
subset(datrep,id==93003 & month==0)[,c('art','cd4_v','threshold','min.allowed','max.allowed')]
subset(datrep,id==93006 & month==0)[,c('art','cd4_v', 'threshold','min.allowed','max.allowed')]

###############################
# Here: Exercise Question 5-6 #
###############################

# 5)
summary(factor(subset(datrep,month==0)$threshold))

# 6) 6 and 4 strategies respectively (makes sense)
subset(datrep,id==93003 & month==0)[,c('art','cd4_v','threshold','min.allowed','max.allowed')]
subset(datrep,id==93006 & month==0)[,c('art','cd4_v', 'threshold','min.allowed','max.allowed')]

#############################
# Let's continue analysis   #
#############################

# CENOSRING: when one clone is deviating from a particular strategy
datrep$c <- rep(0,nrow(datrep))             # indicator for artificial censoring
datrep$elig.c <- rep(NA,nrow(datrep))       # eligibility for artificial censoring

datrep$elig.c <- ifelse(datrep$art_base==1,0,1) # if treatment started as baseline
# now censor (c=1) if ART started before crossing the threshold or
#                  if ART is not started when crossing the threshold
datrep$c <- ifelse(datrep$art==0 & datrep$threshold>datrep$cd4_v & datrep$elig.c==1,1,datrep$c)
datrep$c <- ifelse(datrep$art==1 & datrep$threshold<=datrep$min_cd4_v & datrep$elig.c==1,1,datrep$c)
# Once censored, stay censored
cumsumf <- function(myv){
myv[cumsum(myv)>0] <- 1
return(myv)
}
cumsumf2 <- function(myv){return(cumsum(myv))}
newc <- unsplit(lapply(split(datrep$c, datrep$newid),cumsumf),datrep$newid)
table(newc)
table(ave(datrep$c,datrep$newid,FUN=function(x) ifelse(cumsum(x)>0,1,0)))

datrep$c <- ave(datrep$c,datrep$newid,FUN=function(x) ifelse(cumsum(x)>0,1,0))
datrep$newc2 <- unsplit(lapply(split(datrep$c, datrep$newid),cumsumf2),datrep$newid) 

# Set outcome and censoring indicators missing after artifical censoring (and drop respective lines)
datrep$aidsdeath <- ifelse(datrep$c==1,NA,datrep$aidsdeath)
datrep$censor <- ifelse(datrep$c==1,NA,datrep$censor)
datrep <- subset(datrep,c==0)

# What data remains?
summary(factor(datrep$threshold))
table(datrep$regime,datrep$aidsdeath)

# Examples patients (again):
summary(factor(subset(datrep,id==93003)$threshold))
summary(factor(subset(datrep,id==93006)$threshold))

###############################
# Here: Exercise Question 7-9 #
###############################

# 7) 138
table(datrep$regime,datrep$aidsdeath)

# 8) + 9)
summary(factor(subset(datrep,id==93003)$threshold))
summary(factor(subset(datrep,id==93006)$threshold))

#############################
# Let's continue analysis   #
#############################

# IPTW (see explanation in text)
datrep$denA.line <- ifelse(datrep$eligible2==0,1,	ifelse(datrep$art==0,1-datrep$probA.d, datrep$probA.d))
datrep$dencum <- ave(datrep$denA.line,datrep$newid,	FUN=function(x) cumprod(x))

datrep$w <- 1/datrep$dencum
datrep$w <- ifelse(is.na(datrep$aidsdeath),NA,datrep$w)
datrep <- subset(datrep,is.na(w)==FALSE)

# Truncation
summary(datrep$w)

trunc.cutoff <- quantile(datrep$w,0.99,na.rm=TRUE)
datrep$w.trunc <- ifelse(datrep$w<trunc.cutoff,
	datrep$w,trunc.cutoff)
summary(datrep$w.trunc)

# final MSM
modw.temp <- glm(aidsdeath ~ relevel(factor(regime),'7')
+ month + monthsq
	+ factor(age_0_cat) + SEX + factor(origin) + factor(mode)
	+ factor(year_0_cat) + factor(cd4_0_cat) + factor(rna_0_cat),
	family = binomial(), weights = w.trunc, data = datrep)
summary(modw.temp)
exp(coef(modw.temp)[2:7])

#################################
# Here: Exercise Question 10-12 #
#################################

# ....with valid CI
install.packages('survey')
require(survey)
modw <- svyglm(aidsdeath ~ relevel(factor(regime),'7')
+ month + monthsq
	+ factor(age_0_cat) + SEX + factor(origin) + factor(mode)
	+ factor(year_0_cat) + factor(cd4_0_cat) + factor(rna_0_cat),
	family = quasibinomial(),
	design = svydesign(id = ~id, weights = ~w.trunc, data = datrep))
summary(modw)
exp(coef(modw)[2:7])
exp(confint(modw)[2:7,])
