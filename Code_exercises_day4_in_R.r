#EXERCISE 2 ERASMUS SP 2019

dat_ex2 <- read.csv('/Volumes/DDPS/AdditionalMaterial_R/dat/lab2.csv')

# Check data format
summary(dat_ex2)

View(dat_ex2)

nrow(dat_ex2)

length( unique ( dat_ex2$id ) )

# Example: a few selected patients

subset( dat_ex2 , id == 95001 ) [ , c ('id' , 'month', 'art' , 'censor' , 'death' ) ]
subset( dat_ex2,id == 95002 ) [ , c ('id' , 'month', 'art' , 'censor' , 'death' ) ]
subset( dat_ex2,id == 96119 ) [ , c ('id' , 'month', 'art' , 'censor' , 'death' ) ]

###############################
# Here: Exercise Question 1-4 #
###############################

# 1)

nrow( dat_ex2 )

length( unique ( dat_ex2$id ) )

# 2) at month 9 and 60 patient months:

subset( dat_ex2 , id == 95001 ) [ , c ('id' , 'month', 'art' , 'censor' , 'death' ) ]

# 3) Similarly, read:

subset( dat_ex2,id == 95002 ) [ , c ('id' , 'month', 'art' , 'censor' , 'death' ) ]

subset( dat_ex2,id == 96119 ) [ , c ('id' , 'month', 'art' , 'censor' , 'death' ) ]

#############################
# Let's continue analysis   #
#############################

# Example: pooled logistic as an approximation to Cox model

mod1 <- glm( death ~ art +
                     month +
                     monthsq,
             family = binomial() ,
             data = dat_ex2 )

summary(mod1)

exp( coef( mod1 ) [ 'art' ] )

# adjusted for baseline confounders
mod2 <- glm(death ~ art + month + monthsq
	+ age_0 + I(age_0^2) + SEX + factor(origin) + factor(mode)
	+ year_0 + cd4_0 + I(cd4_0^2) + rna_0 + I(rna_0^2),
	family = binomial(), data = dat_ex2)
summary(mod2)
exp(coef(mod2)['art'])

# adjusted for time-updated covariates
mod3 <- glm(death ~ art + month + monthsq
	+ age_0 + I(age_0^2) + SEX + factor(origin) + factor(mode)
	+ year_0 + cd4_0 + I(cd4_0^2) + rna_0 + I(rna_0^2)
	+ cd4_v + I(cd4_v^2) + rna_v + I(rna_v^2) + aids,
	family = binomial(), data = dat_ex2)
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
	family = binomial(), data = subset(dat_ex2,pastart==0))
	
# nominator  (treatment)
mod2 <- glm(art ~ month + monthsq
	+ age_0 + I(age_0^2) + SEX + factor(origin) + factor(mode)
	+ year_0 + cd4_0 + I(cd4_0^2) + rna_0 + I(rna_0^2),
	family = binomial(), data = subset(dat_ex2,pastart==0))	

#
dat_ex2$probA.d <- ifelse(dat_ex2$pastart==1,1,	predict(mod, type = 'response'))
dat_ex2$probA.n <- ifelse(dat_ex2$pastart==1,1, predict(mod2, type = 'response'))

# probability of NOT being censored
dat_ex2$notcensor <- 1- dat_ex2$censor

# denominator model (censoring)
mod3 <- glm(notcensor ~ art + month + monthsq
	+ age_0 + I(age_0^2) + SEX + factor(origin) + factor(mode)
	+ year_0 + cd4_0 + I(cd4_0^2) + rna_0 + I(rna_0^2)
	+ cd4_v + I(cd4_v^2) + rna_v + I(rna_v^2) + aids,
	family = binomial(), data = dat_ex2)
dat_ex2$probC.d <- predict(mod3, type = 'response')

# nominator model (censoring)
mod4 <- glm(notcensor ~ art + month + monthsq
	+ age_0 + I(age_0^2) + SEX + factor(origin) + factor(mode)
	+ year_0 + cd4_0 + I(cd4_0^2) + rna_0 + I(rna_0^2),
	family = binomial(), data = dat_ex2)
dat_ex2$probC.n <- predict(mod4, type = 'response')

# Calculate stabilized and non-stabilized weights:

# treatment
dat_ex2$A.num <- ifelse(dat_ex2$art==1,dat_ex2$probA.n,1-dat_ex2$probA.n)
dat_ex2$A.den <- ifelse(dat_ex2$art==1,dat_ex2$probA.d,1-dat_ex2$probA.d)
dat_ex2$A.numcum <- ave(dat_ex2$A.num,dat_ex2$id,	FUN=function(x) cumprod(x))
dat_ex2$A.dencum <- ave(dat_ex2$A.den,dat_ex2$id, FUN=function(x) cumprod(x))

dat_ex2$swA <- dat_ex2$A.numcum/dat_ex2$A.dencum
dat_ex2$wA <- 1/dat_ex2$A.dencum

# censoring
dat_ex2$C.numcum <- ave(dat_ex2$probC.n,dat_ex2$id, FUN=function(x) cumprod(x))
dat_ex2$C.dencum <- ave(dat_ex2$probC.d,dat_ex2$id, FUN=function(x) cumprod(x))
dat_ex2$swC <- dat_ex2$C.numcum/dat_ex2$C.dencum
dat_ex2$wC <- 1/dat_ex2$C.dencum

# combined
dat_ex2$sw <- dat_ex2$swA*dat_ex2$swC
dat_ex2$w <- dat_ex2$wA*dat_ex2$wC

# summary
summary(dat_ex2[,c('swA','wA','swC','wC','sw','w')])
dat_ex2$sw.trunc <- ifelse(dat_ex2$sw>10,10,dat_ex2$sw)

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
	family = quasibinomial(), data = dat_ex2, weights = sw.trunc)
summary(modw.temp)
exp(coef(modw.temp)['art'])

# for robust standard errors either use geeglm or survey:

install.packages('survey')
require(survey)
modw <- svyglm(death ~ art + month + monthsq
	+ age_0 + I(age_0^2) + SEX + factor(origin) + factor(mode)
	+ year_0 + cd4_0 + I(cd4_0^2) + rna_0 + I(rna_0^2),
	family = quasibinomial(),
	design = svydesign(id = ~id, weights = ~sw.trunc, data = dat_ex2))
summary(modw)
exp(coef(modw)['art'])
exp(confint(modw)[2,])

