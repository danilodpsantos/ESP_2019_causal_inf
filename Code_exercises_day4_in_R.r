setwd('U:/ESP2019/ESP48/AdditionalMaterial_R/dat') #insert folder path where data are saved locally
dat <- read.csv('lab2.csv')

# Check data format
summary(dat)
head(dat)
nrow(dat)
nrow(dat[!duplicated(dat[,'id']),])

# Example: a few selected patients
subset(dat,id==95001)[,c('id','month','art','censor','death')]
subset(dat,id==95002)[,c('id','month','art','censor','death')]
subset(dat,id==96119)[,c('id','month','art','censor','death')]


###############################
# Here: Exercise Question 1-4 #
###############################

# 1)
nrow(dat)
nrow(dat[!duplicated(dat[,'id']),])

# 2) at month 9 and 60 patient months:
subset(dat,id==95001)[,c('id','month','art','censor','death')]

# 3) Similarly, read:
subset(dat,id==95002)[,c('id','month','art','censor','death')]
subset(dat,id==96119)[,c('id','month','art','censor','death')]

#############################
# Let's continue analysis   #
#############################

# Example: pooled logistic as an approximation to Cox model
mod1 <- glm(death ~ art + month + monthsq,family = binomial(), data = dat)
summary(mod1)
exp(coef(mod1)['art'])

# adjusted for baseline confounders
mod2 <- glm(death ~ art + month + monthsq
	+ age_0 + I(age_0^2) + SEX + factor(origin) + factor(mode)
	+ year_0 + cd4_0 + I(cd4_0^2) + rna_0 + I(rna_0^2),
	family = binomial(), data = dat)
summary(mod2)
exp(coef(mod2)['art'])

# adjusted for time-updated covariates
mod3 <- glm(death ~ art + month + monthsq
	+ age_0 + I(age_0^2) + SEX + factor(origin) + factor(mode)
	+ year_0 + cd4_0 + I(cd4_0^2) + rna_0 + I(rna_0^2)
	+ cd4_v + I(cd4_v^2) + rna_v + I(rna_v^2) + aids,
	family = binomial(), data = dat)
summary(mod3)
exp(coef(mod3)['art'])

###############################
# Here: Exercise Question 5-7 #
###############################

# 5)-> note that crude OR/HR suggests that "ART is bad for you" (because of time-varying confounding affected by prior treatment)
exp(coef(mod1)['art'])
exp(coef(mod2)['art'])
exp(coef(mod3)['art'])

# 6) CD4 cell count i) influences the decision to start treatment and ii) is also affecting survival and is therefore
#    a confounder. CD4 count is also affected by prior treatment decisions (and that's why we can't use simple regression).

# 7) see key answers


#############################
# Let's continue analysis   #
#############################

# IPTW weights

# denominator (treatment)

mod <- glm(art ~ month + monthsq
	+ age_0 + I(age_0^2) + SEX + factor(origin) + factor(mode)
	+ year_0 + cd4_0 + I(cd4_0^2) + rna_0 + I(rna_0^2)
	+ cd4_v + I(cd4_v^2) + rna_v + I(rna_v^2) + aids,
	family = binomial(), data = subset(dat,pastart==0))
	
# nominator  (treatment)
mod2 <- glm(art ~ month + monthsq
	+ age_0 + I(age_0^2) + SEX + factor(origin) + factor(mode)
	+ year_0 + cd4_0 + I(cd4_0^2) + rna_0 + I(rna_0^2),
	family = binomial(), data = subset(dat,pastart==0))	

#
dat$probA.d <- ifelse(dat$pastart==1,1,	predict(mod, type = 'response'))
dat$probA.n <- ifelse(dat$pastart==1,1, predict(mod2, type = 'response'))

# probability of NOT being censored
dat$notcensor <- 1- dat$censor

# denominator model (censoring)
mod3 <- glm(notcensor ~ art + month + monthsq
	+ age_0 + I(age_0^2) + SEX + factor(origin) + factor(mode)
	+ year_0 + cd4_0 + I(cd4_0^2) + rna_0 + I(rna_0^2)
	+ cd4_v + I(cd4_v^2) + rna_v + I(rna_v^2) + aids,
	family = binomial(), data = dat)
dat$probC.d <- predict(mod3, type = 'response')

# nominator model (censoring)
mod4 <- glm(notcensor ~ art + month + monthsq
	+ age_0 + I(age_0^2) + SEX + factor(origin) + factor(mode)
	+ year_0 + cd4_0 + I(cd4_0^2) + rna_0 + I(rna_0^2),
	family = binomial(), data = dat)
dat$probC.n <- predict(mod4, type = 'response')

# Calculate stabilized and non-stabilized weights:

# treatment
dat$A.num <- ifelse(dat$art==1,dat$probA.n,1-dat$probA.n)
dat$A.den <- ifelse(dat$art==1,dat$probA.d,1-dat$probA.d)
dat$A.numcum <- ave(dat$A.num,dat$id,	FUN=function(x) cumprod(x))
dat$A.dencum <- ave(dat$A.den,dat$id, FUN=function(x) cumprod(x))

dat$swA <- dat$A.numcum/dat$A.dencum
dat$wA <- 1/dat$A.dencum

# censoring
dat$C.numcum <- ave(dat$probC.n,dat$id, FUN=function(x) cumprod(x))
dat$C.dencum <- ave(dat$probC.d,dat$id, FUN=function(x) cumprod(x))
dat$swC <- dat$C.numcum/dat$C.dencum
dat$wC <- 1/dat$C.dencum

# combined
dat$sw <- dat$swA*dat$swC
dat$w <- dat$wA*dat$wC

# summary
summary(dat[,c('swA','wA','swC','wC','sw','w')])
dat$sw.trunc <- ifelse(dat$sw>10,10,dat$sw)

#############################
# Here: Exercise Question 8 #
#############################

# 8) Briefly: positivity violations -> small probabilies -> large weights


#############################
# Let's continue analysis   #
#############################

# MSM  (conditional on baseline confounders)

modw.temp <- glm(death ~ art + month + monthsq
	+ age_0 + I(age_0^2) + SEX + factor(origin) + factor(mode)
	+ year_0 + cd4_0 + I(cd4_0^2) + rna_0 + I(rna_0^2),
	family = quasibinomial(), data = dat, weights = sw.trunc)
summary(modw.temp)
exp(coef(modw.temp)['art'])

# for robust standard errors either use geeglm or survey:

install.packages('survey')
require(survey)
modw <- svyglm(death ~ art + month + monthsq
	+ age_0 + I(age_0^2) + SEX + factor(origin) + factor(mode)
	+ year_0 + cd4_0 + I(cd4_0^2) + rna_0 + I(rna_0^2),
	family = quasibinomial(),
	design = svydesign(id = ~id, weights = ~sw.trunc, data = dat))
summary(modw)
exp(coef(modw)['art'])
exp(confint(modw)[2,])

