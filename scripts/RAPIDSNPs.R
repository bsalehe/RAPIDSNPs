##' R scripts for RAPIDSNPs
##' This script describes the analyses for investigating the relationship between genetic variants
##' in the form of single nucleotide polymorphisms (SNPs) and platelet response to agonist activation.
##' The aim is to alternatively examine the SNPs associated with platelet responses using a rapid
##' computational pipeline which involves random Forests (RF) used for screening the important SNPs,
##' followed by stepwsie regression and penalised methods using ridge regression and LASSO.
##' The initial number of trees applied in the random forest model is 500.
##' 
##' 
##' There are four platlet responses involved in the study. These are Fibrinogen binding in 
##' response to ADP (adenosine diphosphate)(FA) and collagen related peptide (CRP-XL)(FC), 
##' and Pselectin expression in response to ADP (PA) and CRP (PC).
##' This script is intended to describe PC but has been applied to other platelet responses as well.
##' 
##'
##' The script was run using R 3.1.2 "Pumpkin Helmet" 64-bit

#==================================================================================================
## Step 1: Loading key packages
#=================================================================================================

## Step 1: Loading packages and libraries
install.packages("randomForest")
install.packages("leaps")
install.packages("ridge")
install.packages("glmnet")
require(ridge)
require(randomForest)
require(leaps)
require(glmnet)

#===================================================================================================

##      Step 2: Data acquisition and Preprocessing: Dataset of 462 subjects and 1430 SNPs

#Load the SNP data for fibrinogen binding in response to ADP (PA)
load_snps_pc <- read.csv("bloodomics_pc.csv", header=TRUE, na.string=c(".", " "))

#Remove firstcolumn of rownames
snp_pc_new <- load_snps_pc[,2:ncol(load_snps_pc)]


#Display some few rows and columns
head(snp_fc_new[,c(1,2,1431)])

#==================================================================================================

##      Step 3: Applying random forest algorithm for screening SNPs

#==================================================================================================

##          500 trees (ntree default)

#==================================================================================================

## Use randomForeest to reduce number of SNPs so that n > p using variable importance scores

## This code run random forest models as regression functions with importance parameter set to be TRUE
## and number of trees (ntree) used in the model are 500 trees, with mtry = p/3 used in partitoning the SNPs in the model

snps.model.pc.rf2 <- randomForest(PselectinCRP ~ ., snp_pc_new, importance = T, ntree = 500,
                                  mtry=(ncol(snp_pc_new)-1)/3)
print(snps.model.pc.rf2)


## Plot the random forest object to visualise the error
plot(1:500, snps.model.pc.rf2$mse, col="green", type="l", xlab="Number of trees", 
     ylab="Test MSE", main="MEAN SQUARED ERROR VIS NUMBER OF TREES WHILE CHANGING MTRY(M) FOR PC MODEL")
legend("topright", c("m=p", "m=p/3", "m=sqrt(p)", "m=25"), col=c("green", "red", "blue", "black"), cex=1, lty=1)


## Plot important snps from snp.model.fc.rf1 RF model
##
## Automatic RF picked SNPs. For this run a total of 30 SNPs was picked
varImpPlot(snps.model.pc.rf1, sort=TRUE, scale=FALSE,
           main=deparse("Important 30 SNPs for PC model"), type=1)

varImpPlot(snps.model.pc.rf1, sort=TRUE, n.var=min(40, nrow(snps.model.pc.rf1$importance)), scale=FALSE,
           main=deparse("Important 40 SNPs for PC model"), type=1)## 40 SNPs were picked


## Obtain and store important snps
imp_snps_rf500_pc <- importance(snps.model.pc.rf1, scale=FALSE)

# Extract top 40 important SNPs approximately sqrt(p)
best40_pc500 <- rownames(imp_snps_rf500_pc)[order(imp_snps_rf500_pc[, "%IncMSE"],
                                                  decreasing = T)[1:40]]
best40_pc500  ## print out selected top 40 SNPs

# New dataset with reduced SNPs embeded with the response variable which is Pselectin in response
# to CRP-XL (pc)
snps40.pc500 <- snp_pc_new[, c(best40_pc500, "PselectinCRP")]

##================================================================================================
## Forward stepwise model
##================================================================================================

# Run the stepwise regression model using the reduced dataset
subs_pc_rf500 <- regsubsets(PselectinCRP ~ ., data = snps40.pc500, method = "forward")
summary(subs_pc_rf500)

#Running linear models to identify the significant SNPs
lmfit_subs_pc_rf500 <- lm(PselectinCRP ~ rs41306982 + rs4504857 + rs41291724 + 
                            rs28514576, data = snps40.pc500)
summary(lmfit_subs_pc_rf500)## Print the linear model
##===============================================================================================
## Ridge regression
##================================================================================================

## Applying ridge regression using the Ridge package
## This is important for finding selected significant SNPs 
## using Wald-test statistics from the ridge regression model
mod_ridge_pc_500 <- linearRidge(PselectinCRP ~ ., data = snps40.pc500)
summary(mod_ridge_pc_500)

##==============================================================================================
## Regression model by lasso methods
#===============================================================================================

## Applying lasso using glmnet package

## Change the dataframe to matrix of SNPs genotypes predictors
geno_pc_rf500 <- as.matrix(snps40.pc500[,1:40])

## Change the response variable which is the Pselectin measure activated by CRP. 
##It is converted into
## a vector of 1 dimension
pheno_pc_rf500 <- snps40.pc500[,41]

## Run the lasso model
lasso_40_pc_rf500 <- glmnet(geno_pc_rf500, pheno_pc_rf500)

## Run cross-validation to identify the best SNPs after shrinking other SNPs
cv_snps40_pc_rf500 <- cv.glmnet(geno_pc_rf500, pheno_pc_rf500, alpha=1, nfolds = 10)
lambda_40_pc_rf500 <- cv_snps40_pc_rf500$lambda.min; # For smallest lambda

res_40_pc_rf500 <- predict(lasso_40_pc_rf500, s=lambda_40_pc_rf500, type="coefficients")
print("LASSO summary model for 40 SNPs after running RF with 500 trees")
## Print all the SNPs showing those with penalised coefficients to zero and others 
##picked by the model.
res_40_pc_rf500 


##' Run the linear model to find the significant SNPs based on the F-test (Partial F-Test)
##' In this case you will be stepwisely selecting those SNPs with fairly large coefficients to
##' the fitted linear model. Selection will be based by looking at 
##' the largest coefficients downward to lower
##' coefficients in the 'res_40_pc_rf500' output. The model might be improved by gradually 
##' removing from the model the SNPs with insignificant p-values ('stepwise variant')
##' This run provided this model
##' 
lmfit_snps_lasso_pc_500 <- lm(PselectinCRP ~ rs41306982 + rs4504857 + rs41291724, snps40.pc500)
summary(lmfit_snps_lasso_pc_500) ## Print out the selected significant SNPs

## =================================================================================================

## Run Boruta methods for further independent SNP selection
library(Boruta)

## The boruta model using selected 40 SNPs
Bor.model.snps40.pc <- Boruta(PselectinCRP ~ ., data=snps40.pc500)
## Print the model
print(Bor.model.snps40.pc)

## A RF wrapper model with Boruta model
model.rf_boruta.snps40.pC <- randomForest(snps40.pc500[,getSelectedAttributes(Bor.model.snps40.pc)],
                                          snps40.pc500$PselectinCRP, importance=T)
print(model.rf_boruta.snps40.pC)

##Plot the selected important snps using boruta plot
par(mar=c(7.2,4.1,2.0,2.1))
plot(Bor.model.snps40.pc, colCode=c("dark gray", "light gray","white","gray"), whichShadow=c(FALSE,FALSE,FALSE,FALSE),
     las = 3, xlab = "", main = "Confirmed SNPs for PC Model")


##Select important and tentative snps
Bor_snps_pc_imp <- getSelectedAttributes(Bor.model.snps40.pc, withTentative = F)
print(Bor_snps_pc_imp)
Bor_snps_pc <- getSelectedAttributes(Bor.model.snps40.pc, withTentative = T)
print(Bor_snps_pc)

# Repeat step 3 with iterations:- 
# 1- ntree = 1000 trees
# 2- ntree = 2000 trees
# 3- ntree = 3000 trees

#============================================================================================================================

##      Step 4: Select the key significant SNPs using consensus approach when the model performance is high (ntree = 3000)

#============================================================================================================================

## Compare all results from the different methods to identify the most significant SNPs by consensus
# summary(lmfit_subs_pc_rf500) # stepwise model (Wald-test)
# summary(mod_ridge_pc_500) # Ridge regression (Wald-Test)
# summary(lmfit_snps_lasso_pc_500) # Lasso using Partial F-test
# cat(Bor_snps_pc_imp) # Boruta method
## Select the consnensu SNPs (Key SNPs)
