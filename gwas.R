# Tutorial https://pbreheny.github.io/adv-gwas-tutorial/quality_control.html 

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

# The .bim file, by contrast, contains information on the genetic loci (SNPs):
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
# you need, in the same directory where you’ve stored the penncath.bed 
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
###### ERRROR! TRY CALL BELOW
system("/Users/neuropromotion/desktop/gwas/plink --bfile /Users/neuropromotion/desktop/gwas/data/penncath --make-bed --out /Users/neuropromotion/desktop/gwas/data/penncath_clean")

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
first_10_snps <- penncath$genotypes[,1:10] # NA values should be filled by imputation:

# Импутация - важный этап обработки генетической информации после QC. ОЧень часто в данных отсутсвует
# часть информации об SNP в виде NA значений. Что бы их заполнить выполняется импутация.
# Pruning - прунинг это метод, который удаляет высоко коррелированные SNPs на основе порога LD (linkage disequilibrium)
# Цель заключается в том, чтобы оставить только независимые (некоррелированные) генетические маркеры для анализа.
# Clumping - это процесс, при котором SNPs группируются (или "кластеризуются") на основе LD с наиболее значимым (по p-value) SNP в каждой области генома.
# Выбирается "ведущий" SNP (обычно с наименьшим p-value) в каждом участке, а затем SNPs, которые находятся в LD с этим "ведущим", кластеризуются и исключаются из дальнейшего анализа.
# 
# Есть статья почему клампинг является лучшей альтернативой прунингу:
# https://privefl.github.io/bigsnpr/articles/pruning-vs-clumping.html
#

######## IMPUTATION
# Tutorial link: https://pbreheny.github.io/adv-gwas-tutorial/imputation.html
# SNPs that are not missing are described as “called.” The call rate is the proportion of genotypes that are called. 
# Therefore, a call rate of 1 indicates that a SNP has no missing values.
# давайте посмотрим на статистику по SNP
obj <- snp_attach(paste0(path,"qc_penncath.rds"))
snp_stats <- bigstatsr::big_counts(obj$genotypes)
colnames(snp_stats) <- obj$map$marker.ID

snp_stats[,1:5]
hist(snp_stats) 

(any_missing <- sum(snp_stats[4,] != 0)) # shows # of SNPs with NA values 
call_rates <- colSums(snp_stats[1:3,])/colSums(snp_stats) 
head(call_rates)

# This tells us that 565138 out of 696644 SNPs still have some missing values (albeit less than 10%.) The average call rate is 0.988 – this high call rate is a good sign, affirming the quality of our data.
# impute based on mode
obj$geno_imputed <- snp_fastImputeSimple(Gna = obj$genotypes,
                                         method = "mode",
                                         ncores = nb_cores())
# save imputed values 
obj$geno_imputed$code256 <- bigsnpr::CODE_IMPUTE_PRED
# look to see that imputation was successful:
imp_snp_stats <- big_counts(X.code = obj$geno_imputed)
imp_snp_stats[,1:5] # all 0s in NA row 

# Let’s save this fully imputed data set for future use in downstream analyses:
obj <- bigsnpr::snp_save(obj)

 
############################################ 
############################################ SUMMARY STATISCTICS
############################################ 
system("/Users/neuropromotion/desktop/gwas/plink --bfile /Users/neuropromotion/desktop/gwas/data/qc_penncath --assoc --out /Users/neuropromotion/desktop/gwas/data/analysis")
# команда для получения файла .log:
# Этот файл содержит журнал выполнения, который включает в себя все сообщения об ошибках, предупреждения и информацию о ходе выполнения команды. Он полезен для диагностики, если возникли проблемы, и для отслеживания настроек, использованных в анализе.
# и файла .assoc:
# Этот файл содержит результаты ассоциативного анализа. В нем будут строки с информацией о SNP, включая: Имя SNP, Хромосому, Позицию, Частоту аллеля,  p-value и другие статистики.
# по сути это SUMMARY STATISTICS GWAS
# Загрузка результатов ассоциативного анализа
assoc_results <- read.table(paste0(path, 'analysis.assoc'), header = TRUE)
assoc_results[1:5,]
assoc_results[assoc_results$P < 1e-5, ]

############################################ HOME ASSIGMENT (22 Chr. Mannh. plot)
## Tutorial: https://pbreheny.github.io/adv-gwas-tutorial/marginal_analysis.html

# reload QC'd & imputed data
obj <- snp_attach(paste0(path, "qc_penncath.rds"))
clinical[1:5,]
# Principal component analysis (PCA)
# SVD 
svd_X <- big_SVD(obj$geno_imputed, # must have imputed data
                 big_scale(), # centers and scales data -- REALLY IMPORTANT! 
                 k = 10 # use 10 PCs for now -- can ask for more if needed
)
# Note: this can take a while to run (several minutes, depending on your
# machine). I'm going to save this SVD as an RDS object, and refer back to that in future use. 

# saveRDS(svd_X, paste0(path, "svd_X.rds"))

# calculate the PCs                          
pc <- sweep(svd_X$u, 2, svd_X$d, "*")
dim(pc) # will have same number of rows as U
# [1] 1401   10
names(pc) <- paste0("PC", 1:ncol(pc))
# Fit a logistic model between the phenotype and each SNP separately
# while adding PCs as covariates to each model
obj.gwas <- big_univLogReg(X = obj$geno_imputed,
                           y01.train = obj$fam$CAD,
                           covar.train = pc[,1:5],
                           ncores = nb_cores())

# this takes a min or two, so I will save these results 
# saveRDS(object = obj.gwas, file = paste0(path, "gwas.rds"))

# Q-Q plot of the object
qq <- bigsnpr::snp_qq(obj.gwas)
(qq)

viridis22 <- c(rep(c("#fde725", "#90d743", "#35b779", "#21918c", "#31688e", 
                     "#443983", "#440154"), 3), "#fde725")

manh <- snp_manhattan(gwas = obj.gwas,
                      infos.chr = obj$map$chromosome,
                      infos.pos = obj$map$physical.pos,
                      colors = viridis22)

(manh)

# NB:  5 × 10e−8 is a common threshold for significance in GWAS studies, 
#   whereas 5 x 10e-6 is a common threshold for "suggestive" results

signif_threshold <- 5e-8 
suggest_threshold <- 5e-6 

(manh + 
    geom_hline(yintercept = -log10(signif_threshold),
               # note: plot y-axis is on -log10 scale 
               color = "#35b779",
               lty = 2) + 
    geom_hline(yintercept = -log10(suggest_threshold),
               color = "#443983",
               lty = 2)
)




# Отфильтровать данные для 22-й хромосомы
chr22_gwas <- obj.gwas[obj$map$chromosome == 22, ]
chr22_pos <- obj$map$physical.pos[obj$map$chromosome == 22]
# Построить манхэттенский график только для 22-й хромосомы
manh <- snp_manhattan(gwas = chr22_gwas,
                      infos.chr = rep(22, length(chr22_pos)), # Только 22-я хромосома
                      infos.pos = chr22_pos,
                      colors = viridis22)
# Отобразить график
(manh + 
    geom_hline(yintercept = -log10(signif_threshold),
               # note: plot y-axis is on -log10 scale 
               color = "#35b779",
               lty = 2) + 
    geom_hline(yintercept = -log10(suggest_threshold),
               color = "#443983",
               lty = 2)
)
# ID ADD and chr info:
obj.gwas$ID <- obj$map$marker.ID
obj.gwas$CHR <- obj$map$chromosome
# Расчет z-статистики
obj.gwas$z_value <- obj.gwas$estim / obj.gwas$std.err
# Преобразование z-статистики в p-value
obj.gwas$p_value <- 2 * pnorm(-abs(obj.gwas$z_value))
# получаем значимые снипы
significant_snps <- obj.gwas[obj.gwas$p_value < 1e-8, ]

genes <- c('Unknown', 'PPP1R12B', 'Unknown', 'CDKN2B-AS1', 'CDKN2B-AS1', 'CDKN2B-AS1',
           'CDKN2B-AS1', 'CDKN2B-AS1', 'CDKN2B-AS1', 'CDKN2B-AS1', 'CDKN2B-AS1')

rownames(significant_snps) <- 1:nrow(significant_snps) # обновим индексы для красоты
significant_snps$mapped_genes <- genes # добавим гены
significant_snps <- significant_snps[, c(5:9)] # почистим 
significant_snps


# GWAS Catalogue
# SS - summary statistics - результат GWAS исследования 
# ID | Chr | Pos | MAF | Major allel | Minor | Effect | p-value 
# Effect это коэффициент B1 перед 
# y = B0 + B1 * X1 + e (где B1 либо 0 (если вариант не альтернативный), 1 (если гетерозигота), 2 (если гомозигота))

# MAF - частота минорного аллеля 
# EAF - эффективная 
# RAF - рисковая аллельная частота 



# Electronic health record (EHR)
# MVP - million veteran programm
# FUMA - fuctional mapping and annotation of genome wide association studies


