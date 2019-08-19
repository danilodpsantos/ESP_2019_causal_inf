## ESP Causal Inference day 3

## R set up:
# Setting read_libs function:
read_libs <- function(packages){

    # Creates a vector with the required but not installed packages:
  install <- packages[!(packages %in% installed.packages()[, "Package"])]
  
  # if there is any uninstalled package in the list:
  
  if(length(install) > 0){
    # Install the uninstalled packages and its dependencies (other required packages in order to function)
    install.packages(pkgs = install, dependencies = TRUE)
  }
  
  # Applying the "require" function in each package of the list, 
  # the "character.only = TRUE" argument is meant to use the packages names as strings
  
    invisible(sapply(packages, require, character.only = TRUE))
}

# Using the read libs function:
# call a vector with the necessary packages names:

read_libs(c("tidyverse", "haven", "survey"))


## Reading data:

dat <- read.csv("/Volumes/DDPS/AdditionalMaterial_R/dat/lab1.csv")

view(dat)

# Look at data
names(dat)

nrow(dat)

view(dat)

# Only use subset of data ("complete cases for outcome")

dat_complete <- filter(dat, 
                       is.na( wt82_71 ) == FALSE)

nrow(dat_complete)

view(dat_complete)

# Descriptives
summary( factor( dat_complete$qsmk ) )

summary( dat_complete$wt82_71 )

mean( subset( dat_complete , qsmk == 1 ) $ wt82_71 )

mean( subset( dat_complete, qsmk==0) $ wt82_71 )

# Create new variable
dat_complete$older <- ifelse( dat_complete$age > 50, 1, 0)

###############################
# Here: Exercise Question 1-4 #
###############################

# 1) Dropped:

nrow(dat) - nrow( dat_complete )

# 2)
mean( dat_complete$wt82_71 )

# 3)
sum( dat_complete$qsmk )

nrow( dat_complete )-sum( dat_complete$qsmk )

# 4)
mean( subset( dat_complete, qsmk==1 )$wt82_71 )

mean( subset( dat_complete, qsmk==0 )$wt82_71 )

#############################
# Let's continue analysis   #
#############################

# weights through 
# 1) tables

table( dat_complete$older,
       dat_complete$qsmk )                   #cell count, columns=qsmk

proportions <- prop.table( table ( dat_complete$older,
                                   dat_complete$qsmk ) , 1 )     #denominator, columns=qsmk

proportions

dat_weigths <- 1/proportions #weights, columns=qsmk

dat_weigths

# 2) logistic regression

mod <- glm( qsmk ~ older, 	family = binomial(), data = dat_complete )

dat_complete$probA.d <- predict( mod, type = 'response' )

dat_complete$A.den <- ifelse(dat_complete$qsmk==1,
                             
dat_complete$probA.d,1-dat_complete$probA.d)

dat_complete$w <- 1/dat_complete$A.den

summary(factor(dat_complete$w))    # ca. 1/prop.table(table(dat_complete$older,dat_complete$qsmk),1) 

# E(Y^a)with weighted mean

weighted.mean( subset( dat_complete, qsmk==1)$wt82_71, w=subset(dat_complete, qsmk==1)$w)

weighted.mean(subset( dat_complete, qsmk==0)$wt82_71, w=subset(dat_complete, qsmk==0)$w)


###############################
# Here: Exercise Question 5-7 #
###############################

# 5)

unique(dat_complete$w)

# 6)
summary(factor( dat_complete$w )) 

dat_weigths

# 7) 
weighted.mean(subset ( dat_complete , qsmk == 1 )$wt82_71,
              w = subset ( dat_complete , qsmk == 1 ) $ w ) - 
  
weighted.mean(subset ( dat_complete , qsmk == 0 ) $ wt82_71,
              w=subset(dat_complete,qsmk==0)$w )


#############################
# Let's continue analysis   #
#############################

# IPTW - standard weights

A.den <- glm ( qsmk ~ sex + 
                      factor( race ) +
                      age +
                      I( age^2 ) +
                      factor( education ) +
                      smokeintensity +
                      I( smokeintensity^2 ) +
                      smokeyrs +
                      I( smokeyrs^2) +
                      factor(exercise) +
                      factor(active) +
                      wt71 +
                      I(wt71^2), 
                family = binomial() ,
                data = dat_complete)
	
dat_complete$probA.den <- predict( A.den, type = "response" )

dat_complete$w <- ifelse(dat_complete$qsmk == 1,
                         1/ dat_complete$probA.den,
                         1/( 1 - dat_complete$probA.den ) )

# IPTW - stabilized weights	

A.num <- glm( qsmk ~ 1, 
              family = binomial(),
              data = dat_complete )

dat_complete$probA.num <- predict( A.num, type = "response" )

dat_complete$sw <- ifelse(dat_complete$qsmk == 1,
                          dat_complete$probA.num / dat_complete$probA.den,
                          ( 1 - dat_complete$probA.num ) / ( 1 - dat_complete$probA.den ) )

# Summary of weights

summary( dat_complete [ , c('w','sw') ] )


###############################
# Here: Exercise Question 5-7 #
###############################

# 8) + 9)

summary( dat_complete [ , c('w','sw') ] )

#############################
# Let's continue analysis   #
#############################

#  mean outcomes in the weighted population

weighted.mean( subset( dat_complete , qsmk == 1 )$wt82_71 , w = subset( dat_complete, qsmk == 1 )$w )
weighted.mean( subset( dat_complete , qsmk == 0 )$wt82_71 , w = subset( dat_complete, qsmk == 0 )$w)

weighted.mean( subset( dat_complete, qsmk == 1 )$wt82_71 , w = subset( dat_complete , qsmk == 1 )$sw)
weighted.mean( subset( dat_complete, qsmk == 0 )$wt82_71 , w = subset( dat_complete , qsmk == 0 )$sw)


#############################
# Let's continue analysis   #
#############################

# final MSM analysis
modw <- glm( wt82_71 ~ qsmk ,
             data = dat_complete ,
             weights = sw)

summary(modw)

# for better confidence intervalsuse package "survey"

modw <- svyglm(wt82_71 ~ qsmk,
               design = svydesign ( id = ~ seqn,
                                    weights = ~ sw,
                                    data = dat_complete ) )

summary(modw)

confint(modw)


#########################################
# now the g-formula (standardization)   #
#########################################

# prepare dataset
dat1 <- dat_complete #first copy of dataset (equal to original)

dat1$interv <- rep(-1 , nrow( dat1 ) ) 

dat2 <- dat_complete #second copy of dataset (set tx to 0, outcome to missing)

dat2$interv <- rep( 0 , nrow( dat2 ) )

dat2$qsmk <- rep( 0 , nrow( dat2 ) )

dat2$wt82_71 <- rep( NA , nrow( dat2 ) )

dat3 <- dat_complete #third copy of dataset (set tx to 1, outcome to missing)

dat3$interv <- rep( 1 , nrow( dat3 ) )

dat3$qsmk <- rep( 1 , nrow( dat3 ) )

dat3$wt82_71 <- rep( NA , nrow( dat3 ) )

onesample <- as.data.frame( rbind( dat1 , dat2 , dat3 ) )

dim( onesample )

summary( factor( onesample$interv ) )

# Preparing the linear mmodel

modst <- glm(wt82_71 ~ qsmk +
                       sex +
                       factor(race)	+
                       age +
                       I(age^2) +
                       factor(education) +
                       smokeintensity +
                       I(smokeintensity^2) +
                       smokeyrs +
                       I(smokeyrs^2) +
                       factor(exercise) +
                       factor(active) +
                       wt71 +
                       I(wt71^2), 
             data = onesample)

summary(modst)

# Couldmzt understand what did they do here :)
onesample$meanY <- predict( modst , onesample , type = "response" )

with(onesample,
     tapply( meanY , list( interv ), mean ) )





