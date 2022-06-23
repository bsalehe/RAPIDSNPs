## This script is used for reproducing the data for validating the pipeline RAPIDSNPs.
## The data contains artificially generated 1400 SNPs with 460 randomly sampled subjects.
## The phenotype or response variable is normally distrubuted with mean 0 and 1 standard deviation.
#

## Set seed for reproducing the data
#
set.seed(12345)

## Simulating the artificial SNPs with their genotypes (independent variables) 1400 SNPs, and 460 subjects
## Creating the data frame of 460 subjects and 1400 SNPs (df_SNPs)
#
df_SNPs <- data.frame(replicate(1400, sample(1:3, 460, rep=TRUE)))

## Simulating response variable (phenotype) with random normal distribution.
#
yres_phen <- rnorm(460)


## Combining the response variable and SNPs, i.e, adding a new column to the above created data frame of SNPs (df_SNPs)
#
df_SNPs["phenotype"] <- yres_phen
df_SNPs_yres_phen <- df_SNPs[, 1:ncol(df_SNPs)] ## New data set with the SNPs genotypes and phenotype.


## Viewing the few first rows and columns of the created data frame
#
head(df_SNPs_yres_phen[,c(1:4,1401)])


