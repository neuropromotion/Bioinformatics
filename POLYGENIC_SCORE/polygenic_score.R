library(data.table)
library(dplyr)
library(ggplot2)
library(bigsnpr)
library(pROC)

path <- "YOUR PATH TO DATA" 
# данные по фенотипу
clinical <- fread(paste0(path, 'penncath.csv'))
head(clinical)
# прошедшие QC и импутации данные по генетике:
penncath <- snp_attach(paste0(path, 'qc_penncath.rds')) 

# файл с снипами и их beta
file_path <- "PATH TO PGS FILE" 

# Читаем файл, пренебрегая всеми комментами
pgs.df <- read.table(file_path, header = F, comment.char = "#")
pgs.df[1:5, ] # смотрим как получилось
pgs.df <- pgs.df[, c(1, 4, 6)] # сохраняем только то что нужно

head(penncath$map) # смотрим какие столбцы нас интересуют (это 2 - ID снипа, 5 и 6 - нуклеотиды в аллелях)
aggr.df <- penncath$map[, c(2,5:6)] # отбираем только их
head(aggr.df) # смотрим как получилось

table(clinical$CAD) # выходит 468 случаев без CAD и 933 с CAD
clinical$ind <- 1:nrow(clinical) # добавляем индексы как новую колонку (нужно дальше)

################################################
################################################
# теперь надо найти все снипы которые пересекаются у наших пациентов и в файле pgs:
indicies <- vector() # создаем вектор куда будем записывать соответствующие индексы
SNP.intersection <- vector() # и запомним все снипы которые есть в пересечении

# идем циклом
# ctr <- 0
#for (i in 1:length(pgs.df$V1)){
#  ctr <- ctr + 1
#  id <- pgs.df$V1[i] # запоминаем текущий ID снипа 
#  ind <- which(aggr.df$marker.ID == id) # запоминаем индекс снипа с которых он смэтчился
#  # если данные по такому снипу есть у нас, то length(ind) != 0, запоминаем индекс в вектор
#  if (length(ind) != 0){
#    indicies <- c(indicies, ind)
#    SNP.intersection <- c(SNP.intersection, id)
#    ind <- NULL # обновляем значение ind
#  }
#  if (ctr %% 10000 == 0){
#    print(as.character(ctr), '/', as.character(length(pgs.df$V1)))
#  }
#}
################################ можно в миллион раз проще, это же R ########################
################################
result.df <- pgs.df[pgs.df$V1 %in% aggr.df$marker.ID, ]
result.df[1:5,]
indicies.pgs.df <- which(pgs.df$V1 %in% aggr.df$marker.ID) # получаем все индексы снипов из пересечения в pgs.df
# indicies <- match(SNP.intersection, pgs.df$V1) # второй способ
SNP.intersection <- result.df$V1 # получаем все ID снипов в пересечении
indicies.aggr.df <- which(aggr.df$marker.ID %in% SNP.intersection) # получаем все индексы снипов из пересечения в aggr.df

# assert
if (length(indicies) != length(SNP.intersection)){
  print('ERROR')
}

################
# запоминаем индексы всех пациентов с CAD и без
zeros.inds <- clinical[clinical$CAD == 0]$ind
ones.inds <- clinical[clinical$CAD == 1]$ind
# создаем две матрицы снипов из наших данных (для пациентов с CAD и без c теми снипами которых пересеклись). 
# penncath$genotypes - матрица снипов, по строкам
# пациенты, по столбцам снипы (в нашем случае 1401х696644). 
# dim(penncath$genotypes)

zeros.mat <- penncath$geno_imputed[zeros.inds, indicies.aggr.df]
ones.mat <- penncath$geno_imputed[ones.inds, indicies.aggr.df]
dim(zeros.mat)# получается матрица 468 х k (для пациентов без CAD)

# теперь создаем две нулевые матрицы для каждого пациента, будем их заполнять скорами
values.for.ones <- matrix(0, nrow = dim(ones.mat)[1])
values.for.zeros <- matrix(0, nrow = dim(zeros.mat)[1])


# основной цикл
# Получаем данные об эффектном аллеле и весах из pgs.df для снипов из пересечения:
alleles_pgs <- pgs.df$V4[indicies.pgs.df]
weights_pgs <- pgs.df$V6[indicies.pgs.df]

for (i in seq_along(SNP.intersection)) { # итерируемся по снипам из пересечения
  id <- SNP.intersection[i] # берем ID снипа
  allel <- alleles_pgs[i] # по этому индексу берем эффектную аллель
  weight <- weights_pgs[i] # и вес 
  
  # теперь берем первую и вторую аллель наших пациентов:
  major <- aggr.df[indicies.aggr.df[i], ]$allele1 
  minor <- aggr.df[indicies.aggr.df[i], ]$allele2 
  # если минорная аллель соответствует эффектной алелли из PGS, то:
  # значит 0 в матрице генотипов значит гомозиготность по минорной аллели
  # меняем нули на двойки, двойки на нули, единицы (гетерозиготность) оставляем
  # домножаем на коэффициент бета и прибавляем к итоговому скору пациента
  if (minor == allel){
    zeros.mat[, i][zeros.mat[, i] == 0] <- -1 # временно меняем нули на -1
    zeros.mat[, i][zeros.mat[, i] == 2] <- 0 # двойки на нули
    zeros.mat[, i][zeros.mat[, i] == -1] <- 2
    values.for.zeros <- values.for.zeros + (zeros.mat[, i] * weight)
    
    ones.mat[, i][ones.mat[, i] == 0] <- -1
    ones.mat[, i][ones.mat[, i] == 2] <- 0
    ones.mat[, i][ones.mat[, i] == -1] <- 2 
    values.for.ones <- values.for.ones + (ones.mat[, i] * weight)
  }
  if (major == allel){ 
    values.for.zeros <- values.for.zeros + (zeros.mat[, i] * weight)
    values.for.ones <- values.for.ones + (ones.mat[, i] * weight)
  }
  
  # Печатаем каждые 10 000 итераций
  if ((i %% 10000) == 0) {
    print(paste0(as.character(i), '/', as.character(length(SNP.intersection))))
  }
}

# нормализуем
norm_ones <- values.for.ones / dim(zeros.mat)[2]
norm_zeros <- values.for.zeros / dim(zeros.mat)[2]

ggplot() +
  geom_density(aes(x = norm_ones[,1], fill = "With CAD"), alpha = 0.5, color = "black") +
  geom_density(aes(x = norm_zeros[,1], fill = "Without CAD"), alpha = 0.5, color = "black") +
  scale_fill_manual(name = "Groups", values = c("With CAD" = "red", "Without CAD" = "green")) +
  theme_minimal() +
  labs(title = "Density plot", x = "Value", y = "Frequency") +
  theme(legend.position = "top")


### Теперь немного преобразуем данные что бы нарисовать ROC-AUC
# Преобразуем матрицы в датафреймы
df.zeros <- data.frame(values = norm_zeros[,1])
df.ones <- data.frame(values = norm_ones[,1])

# Добавляем метку (0 для первой матрицы и 1 для второй)
df.zeros$label <- 0
df.ones$label <- 1

# Объединяем оба датафрейма в один
combined.df <- rbind(df.zeros, df.ones)
head(combined.df)

# Создание ROC-объекта
roc_obj <- roc(combined.df$label, combined.df$values)

# Построение ROC-кривой
plot(roc_obj, 
     main="ROC-кривая", 
     col="blue", 
     lwd=2, 
     print.auc=TRUE, 
     print.auc.x=0.6, 
     print.auc.y=0.2)


# LDlink
# LD clumping - метод для того что бы убрать сцепленные спины
# LD score регрессия
# посчитать генетическую корреляцию между двумя признаками














