library(data.table)
library(dplyr)
library(ggplot2)
library(bigsnpr)

# фильтруем по MAF < 0.05
system('/Users/neuropromotion/desktop/gwas/plink --bfile /Users/neuropromotion/desktop/gwas/clustering/biengi/biengi --maf 0.05 --make-bed --out /Users/neuropromotion/desktop/gwas/clustering/biengi_filtered/biengi_maf_filtered') 
# prunning
system('/Users/neuropromotion/desktop/gwas/plink --bfile /Users/neuropromotion/desktop/gwas/clustering/biengi_filtered/biengi_maf_filtered --indep-pairwise 50 5 0.2 --out /Users/neuropromotion/desktop/gwas/clustering/biengi_prunned/biengi_prunned')
# собираем обратно bed из прунненых файлов
system('/Users/neuropromotion/desktop/gwas/plink --bfile /Users/neuropromotion/desktop/gwas/clustering/biengi_filtered/biengi_maf_filtered --extract /Users/neuropromotion/desktop/gwas/clustering/biengi_prunned/biengi_prunned.prune.in --make-bed --out /Users/neuropromotion/desktop/gwas/clustering/prunned_bed/biengi_pruned')

path <- "/Users/neuropromotion/desktop/gwas/clustering/prunned_bed/" 

bed_file <- snp_readBed(bedfile = paste0(path, 'biengi_pruned.bed'))

data <- snp_attach(paste0(path, 'biengi_pruned.rds')) 

## Проводим импутацию данных:
snp_stats <- bigstatsr::big_counts(data$genotypes)
colnames(snp_stats) <- data$map$marker.ID
snp_stats[,1:5]
hist(snp_stats) 
data$geno_imputed <- snp_fastImputeSimple(Gna = data$genotypes,
                                         method = "mode",
                                         ncores = nb_cores())

## создаем матрицу куда скопируем матрицу фичей:
mat <- matrix(0, nrow = nrow(data$geno_imputed), ncol = ncol(data$geno_imputed))
for (i in 1:dim(data$geno_imputed)[2]){
  mat[,i] <- data$geno_imputed[,i]
}

dim(mat)
## центрируем что бы в последствии перемножить матрицу на себя транспонированную
## и получить матрицу ковариации в соответствии с формулой ковариации:
## cov(X,Y) = 1/(n-1) * ∑((X_i - Xˆ) * (Y_i - Yˆ))
## центрируем как раз для того что бы вычесть из кажого наблюдения - Xˆ - среднюю по признаку 
mat <- scale(mat, center = TRUE, scale = FALSE)
cov <- (mat %*% t(mat)) / dim(mat)[2] # матрица ковариации
# получаем собственные вектора и их значения из матрицы, где:
# собственный вектор это прямая в пространстве куда были спроецированы данные
# а собственное значение это доля обьясненной дисперсии (это свойство мы используем что бы нарисовать далее elbowplot)
eigen_result <- eigen(cov) 

# рисуем elbow, возьмем первые 15 собственных вектором и их значений и посчитаем сколько
# процентов дисперсии они обьясняют
elbow <- vector()
for (i in 1:15){
  elbow <- c(elbow, eigen_result$values[i]/sum(eigen_result$values)*100)
}
elbow_data <- data.frame(
  Component = 1:15,
  Variance = elbow
)
# Рисуем график elbow
ggplot(elbow_data, aes(x = Component, y = Variance)) +
  geom_point(color = "blue", size = 3) +          # Точки
  geom_line(color = "blue", size = 1) +          # Линии
  labs(
    title = "Scree Plot (Elbow Method)",         # Заголовок
    x = "Component Number",                      # Подпись оси X
    y = "Explained Variance (%)"                # Подпись оси Y
  ) +
  theme_minimal() +                              # Минималистичная тема
  theme(
    plot.title = element_text(hjust = 0.5)      # Центрируем заголовок
  )

# Проекция данных на первые две главные компоненты
pc1 <- t(cov) %*% eigen_result$vectors[, 1]
pc2 <- t(cov) %*% eigen_result$vectors[, 2]
pc3 <- t(cov) %*% eigen_result$vectors[, 3]
pc4 <- t(cov) %*% eigen_result$vectors[, 4]
# Объединяем в датафрейм
pca_df <- data.frame(PC1 = pc1, PC2 = pc2, PC3 = pc3, PC4 = pc4)

# можно рисовать график PC1-PC2
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(color = "darkgreen", size = 0.7) +
  labs(
    title = "PCA Plot",
    x = "Principal Component 1",
    y = "Principal Component 2"
  ) +
  theme(
    panel.background = element_rect(fill = "lightgray", color = "black"), # Фон панели
    plot.background = element_rect(fill = "lightgray", color = "black")    # Общий фон
  )

 
# добавляем национальности 

pedind_data <- read.table("/Users/neuropromotion/desktop/gwas/clustering/biengi/biengi.pedind", header = FALSE)
head(pedind_data) # нам нужен pedind_data$Nation

pca_df$Nation <- pedind_data$Nation

ggplot(pca_df, aes(x = PC1, y = PC2, color = Nation)) +
  geom_point(size = 0.7) +
  labs(
    title = "PCA Plot",
    x = "Principal Component 1",
    y = "Principal Component 2"
  ) +
  theme(
    panel.background = element_rect(fill = "lightgray", color = "black"), # Фон панели
    plot.background = element_rect(fill = "lightgray", color = "black")    # Общий фон
  )



##################### Alternative 
pca_result <- prcomp(mat, center = TRUE, scale. = TRUE)
pca_df <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2])

ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(color = "grey", size = 1) +
  labs(
    title = "PCA Plot",
    x = "Principal Component 1",
    y = "Principal Component 2"
  ) +
  theme_minimal()












