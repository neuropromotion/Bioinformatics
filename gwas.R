# Tutorial https://pbreheny.github.io/adv-gwas-tutorial/quality_control.html#Downloading_PLINK

#install.packages('data.table')
#install.packages("bigsnpr")
library(data.table)
library(dplyr)
library(ggplot2)
library(bigsnpr)

path <- "/Users/neuropromotion/desktop/gwas/data/"
path. <- "/Users/neuropromotion/desktop/gwas/"

# fam contains information on the subjects
# The six colums are: Family ID, Individual ID, Paternal ID, Maternal ID, 
# Sex (1=male; 2=female; other=unknown), Phenotype
fam <- fread(paste0(path, 'penncath.fam'))
head(fam)

# Phenotype is typically used to record case-control status or something like that, 
# but it is also quite common to just record clinical/biological information in a separate 
# spreadsheet, which is what was done here:
clinical <- fread(paste0(path, 'penncath.csv'))
head(clinical)

# As you can see, we’ve got the FamID to match this spreadsheet up with the genetic data, 
# the disease status (CAD=1 means that the subject has coronary artery disease), 
# and some covariates (age, triglycerides, HDL and LDL cholesterol levels).

#The .bim file, by contrast, contains information on the genetic loci (SNPs):
bim <- fread(paste0(path, 'penncath.bim'))
head(bim)
# The columns are:
# chromosome (1-22, X, Y or 0 if unplaced)
# rs# or snp identifier
# Genetic distance (morgans)
# Base-pair position (bp units)
# Allele 1 (usually minor allele, or ‘effect’ allele)
# Allele 2 (usually major allele)
############################### 
# So, for example, the file tells us that genetic locus rs12565286 is located 721290 bases 
# into chromosome 1, and that most people have a C there, but some have a G. Of course, 
# “most people” will be relative to a population

bim[bim$V1 == 22] # our chromosome to be analysed

# Finally, the .bed file, which has all the data. 
# To read this into R, run the following chunk only once. This will create the .rds and .bk files
# you need, in the same directory where you’ve stored the penncath.bed (my directory is called data/)
snp_readBed(paste0(path, 'penncath.bed'))

# Now, read in the data to your R session:
penncath <- snp_attach(paste0(path, 'penncath.rds'))

# check out this bigSNP object 
str(penncath)
# Then we have to provide quality control:
# One helpful first step when getting any genetic data is to address the question: 
# what chromosomes are the SNPs on? In our data, we can see this in the map element 
# of our bigSNP object:
penncath$map$chromosome |> table()
# We notice here that (1) we do not have any sex chromosomes2, and 
# (2) we do not have any chromosomes with numbers outside 1-22. 
# Sometimes, SNPs with unknown chromosome are labeled with a chromosome of 0.
# Another common convention is to label mitochondrial variants as having chromosome ‘26’.

############################### Then we have to check Missing data

# Any SNP with a lot of missing data are often excluded from analysis.
# Looking at our genotype data (i.e., what is in our .bed file), 
# we see that some of our SNPs have missing values.
?bigstatsr::big_counts # handy function to summarize data 
big_counts(penncath$genotypes, ind.row = 1:10, ind.col = 1:10) # we see that 4th row has lot of zeros

# let's recieve statistics:
snp_stats <- big_counts(penncath$genotypes)
dim(snp_stats)

# boxplot of 4th row:
boxplot(snp_stats[4,])
summary(snp_stats[4,])# There are definitely some SNPs with lots of missing values – 
# at least one SNP is missing across all 1,401 samples!
# A common practice is to exclude SNPs with >10% missing data. 

# We also need to consider any samples that have high proportions of missing values:
sample_stats <- big_counts(penncath$genotypes, byrow = TRUE)
sample_stats[, 1:10]
boxplot(sample_stats[4,])
summary(sample_stats[4,])

############################### Heterozygosity check
# If an individual had a ton of A/B calls but no A/A or B/B calls, or vice versa, 
# that would likely indicate something was wrong in that sample – we would typically 
# expect all samples to have some heterozygous calls and some homozygous calls.


allele_dat <- sweep(x = sample_stats[1:3,], 
                    # ^ leave off 4th row -- don't count NAs here
                    MARGIN = 1,
                    STATS = c(0, 1, 0),
                    FUN = "*") |> colSums()

boxplot(allele_dat/ncol(snp_stats))
hist(allele_dat/ncol(snp_stats),
     main = "Zygosity",
     xlab = "Proportion of samples which are heterozygous") # should be bell-curved shaped

# No big outliers here… that’s a good sign.

############################### Minor allele frequency (MAF) filtering
# we will need to exclude SNPs that have low variation
# We see there are some SNPs for which all 1401 samples have the A/A alleles, or 
# the major/major alleles. For analysis purposes, we typically need to exclude such SNPs 
# in which variation is rare.

hist(snp_stats[1,])
summary(snp_stats[1,])
# The histogram and summary shows that there are some SNPs for which few, if any, samples have a minor allele. We will need to filter out these rare variants.

############################### ############################### 
############################### Summary ############################### 
# In summary, here are the basic QC steps that (almost) all GWAS need to consider:
# 1) We need to check the chromosomes we have represented in our data

# 2) We need to filter out samples that have a lot of missing values (e.g., samples that had poor quality in upstream data collection/analysis)

# 3) We need to check for heterozygosity, and filter out samples that are clear outliers

# 4) We need to filter out SNPs that have a really low minor allele frequency

# 5) We need to filter out SNPs that are far outside Hardy-Weinberg equilibrium

# 6) We need to consider relatedness in our data – both the relatedness that is known as well as any latent relationship structure(s).
############################### ############################### 
############################### ############################### 





path_to_qc_penncath <- snp_plinkQC(
  plink.path = paste0(path., 'plink'), # you will need to change this according to your machine!
  prefix.in = paste0(path, 'penncath'), # input current data
  prefix.out = paste0(path, 'qc_penncath'), # creates *new* rds with quality-controlled data
  maf = 0.01
)
#system("plink --bfile data/penncath --make-bed --out data/penncath_clean")

path_to_qc_penncath_bed <- snp_plinkQC(
  plink.path = paste0(path., 'plink'), # again, you may need to change this according to your machine!
  prefix.in = paste0(path, 'penncath_clean'), # input data
  prefix.out = paste0(path, 'qc_penncath'), # creates *new* rds with quality-controlled data
  maf = 0.01, # filter out SNPs with MAF < 0.01
  geno = 0.1, # filter out SNPs missing more than 10% of data
  mind = 0.1, # filter out SNPs missing more than 10% of data,
  hwe = 1e-10, # filter out SNPs that have p-vals below this threshold for HWE test
  autosome.only = TRUE # we want chromosomes 1-22 only
  
)
# Let’s look at our new RDS object, the one we just created using our QC steps:

qc_file_path <- snp_readBed(bedfile = paste0(path, 'qc_penncath.bed'))

qc <- snp_attach(qc_file_path)

dim(qc$genotypes)


# Recieving 22 chromosome data
chromosome_22 <- qc$map[qc$map$chromosome == 22, ]

#./plink --bfile /Users/neuropromotion/desktop/gwas/data/qc_penncath --assoc --out /Users/neuropromotion/desktop/gwas/data/new

# команда для получения файла .log:
# Этот файл содержит журнал выполнения, который включает в себя все сообщения об ошибках, предупреждения и информацию о ходе выполнения команды. Он полезен для диагностики, если возникли проблемы, и для отслеживания настроек, использованных в анализе.
# и файла .assoc:
# Этот файл содержит результаты ассоциативного анализа. В нем будут строки с информацией о SNP, включая: Имя SNP, Хромосому, Позицию, Частоту аллеля,  p-value и другие статистики.

# Загрузка результатов ассоциативного анализа
assoc_results <- read.table(paste0(path, 'new.assoc'), header = TRUE)

# Просмотр первых строк результатов
head(assoc_results)
head(assoc_results[assoc_results$CHR == 22, ])

# Извлечение p-values
p_values <- assoc_results$P
p_values_chr22 <- assoc_results[assoc_results$CHR == 22, ]$P

# Добавляем значения p к данным о SNP
qc$map$p <- p_values
chromosome_22$p <- p_values_chr22

# Подготовка данных для манхэттен-графика
chromosome_22$position <- 1:nrow(chromosome_22)

# Построение манхэттен-графикаs
ggplot(chromosome_22, aes(x = position, y = -log10(p))) +
  geom_point(color='lightgreen',alpha = 0.5) +
  theme_minimal() +
  theme_dark() +
  labs(title = "Manhattan Plot for Chromosome 22",
       x = "Position on Chromosome 22",
       y = "-log10(p-value)") +
  theme(axis.text.x = element_blank())


# Подсчет SNP с p-value больше 1e-5
qc$map[qc$map$p < 1e-5, ]











