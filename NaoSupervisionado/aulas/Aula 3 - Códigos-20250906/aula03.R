##################################################################
# Mineracao de Dados Complexos -- MDC 
# Aprendizado de Máquina Não Supervisionado
# Prof. Helio Pedrini
# Codigos da Aula 3 - Reducao de Dimensionalidade
##################################################################




# Fatoração de Matrizes
# Decomposição em Valores Singulares

# cria matriz de entrada
X = as.matrix(data.frame(c(3,2,7,2,9),
                         c(1,9,6,1,7),
                         c(7,8,8,4,6),
                         c(5,4,6,6,3)))

# aplica decomposição em valores singulares
X_svd <- svd(X)

# verifica resultado da multiplicação das três matrizes
X_svd$u %*% diag(X_svd$d) %*% t(X_svd$v)

# U^T * U = I
t(X_svd$u) %*% X_svd$u

# V^T * V = I
t(X_svd$v) %*% X_svd$v

# aproxima matriz de entrada com 3 fatores
k <- 3
X_app <- X_svd$u[,1:k] %*%
         diag(X_svd$d[1:k]) %*%
         t(X_svd$v[,1:k])

# Fatoração de Matrizes Nao-Negativas
# carrega biblioteca
library(nnls)
source('auxiliar_aula3.R')
# cria matriz de entrada
X = as.matrix(data.frame(c(3,2,7,2,9),
                         c(1,9,6,1,7),
                         c(7,8,8,4,6),
                         c(5,4,6,6,3)))

# aplica decomposição NMF com k=3
res <- rnmf(X, 3)

# verifica matrizes W e H
W <- res$W
H <- res$H

# verifica resultado da aproximação da matriz de entrada
W %*% H

# Decomposição CUR de Matrizes
# carrega biblioteca
library(rsvd)

# cria matriz de entrada
X = as.matrix(data.frame(c(3,2,7,2,9),
                         c(1,9,6,1,7),
                         c(7,8,8,4,6),
                         c(5,4,6,6,3)))

# aplica decomposição CUR com k=3
res <- rcur(X, 3)

# verifica matrizes W e H
C <- res$C
U <- res$U
R <- res$R

# verifica resultado da aproximação da matriz de entrada
C %*% U %*% R


# Redução de Dimensionalidade no Pacote R
# Carrega as bibliotecas e o conjunto de dados:
library(mlbench)
library(caret)
data(PimaIndiansDiabetes)

# Apresenta um resumo da base:
summary(PimaIndiansDiabetes)

# Calcula a matriz de correlacao:
matriz_corr <- cor(PimaIndiansDiabetes[,1:8])

# Mostra a matriz de correlacao:
print(matriz_corr)

# Encontra atributos mais relacionados (por exemplo: > 0.50)
altamente_corr <- findCorrelation(matriz_corr, cutoff = 0.5)

# Mostra os indices dos atributos altamente correlacionados:
print(altamente_corr)

# Carrega as bibliotecas e o conjunto de dados:
library(mlbench)
library(caret)
library(randomForest)
data(PimaIndiansDiabetes)

# Para garantir repetibilidade dos resultados
set.seed(7)

# Define o controle usando uma funcao de selecao com florestas aleatorias
controle <- rfeControl(functions = rfFuncs, method = "cv", number =10)

# Executa o algoritmo RFE:
resultados <- rfe(PimaIndiansDiabetes[,1:8], PimaIndiansDiabetes[,9],
                  sizes = c(1:8), rfeControl = controle)

# Mostra os resultados:
print(resultados)

# Apresenta os atributos escolhidos:
predictors(resultados)

# Exibe os resultados:
plot(resultados, type=c("g", "o"))

# Carrega o conjunto de dados:
library(datasets)
data(iris)

# Apresenta um resumo da base:
summary(iris)

# Calcula as componentes principais:
iris.pca1 <- prcomp(iris[,1:4], scale.=TRUE)
iris.pca1

# Sem o uso de scale.=TRUE:
iris.pca2 <- prcomp(iris[,1:4])
iris.pca2

# Apresenta resumo das componentes principais:
summary(iris.pca1)

# Apresenta resumo das componentes principais:
summary(iris.pca2)

# O comando princomp gera os mesmos resultados, mas a chamada e
# um pouco diferente:
iris.pca3 <- princomp(iris[,1:4], cor=TRUE)
iris.pca3
summary(iris.pca3)

# Mostra as primeiras componentes principais:
head(iris.pca1$x)

# O comando iris.pca1$x[,1:2] mostra os dados com apenas 2
# dimensoes nas coordenadas das componentes principais:
head(iris.pca1$x[,1:2])

# Conversao dos dados de volta para as dimensoes originais:
z <- iris.pca1$x[,1:2] %*% t(iris.pca1$rotation[,1:2])
head(z)

# Resumo das componentes principais:
summary(iris.pca1)

# Visualizacao dos componentes principais e respectivas variancias:
screeplot(iris.pca1, type="lines")

# Outra forma de calcular e apresentar a contribuicao de cada
# componente principal:
# calcula desvio padrao de cada componente principal
pr_std_dev <- iris.pca1$sdev

# calcula variancia
pr_var <- pr_std_dev^2

# mostra variancia das 4 componentes principais
pr_var

# Proporcao da variancia explicada:
prop_varex <- pr_var/sum(pr_var)
prop_varex

# Visualiza proporcao da variancia explicada:
plot(prop_varex, xlab="Componentes Principais",
     ylab = "Proporcao de Variancia Explicada", type = "b")

# Visualiza variancia acumulada:
plot(cumsum(prop_varex), xlab = "Componentes Principais",
     ylab = "Variancia Explicada Acumulada", type = "b")

# Exibicao de como os atributos serao transformados:
biplot(iris.pca1, xlab="Componente Principal 1",
       ylab="Componente Principal 2", scale=0)

# Visualizacao de dados com dois atributos:
plot(iris[,1:2], col=iris[,5])

# Visualizacao de componentes principais:
plot(iris.pca1$x[,1:2], col=iris[,5])

# # Leitura arquivo da base wine versão online:
# wine <- read.table(
#   "http://archive.ics.uci.edu/ml/machine-learning-databases/wine/wine.data",
#   sep=",")

# Leitura arquivo da base wine com arquivo offline
# baixado do Moodle:
wine <- read.table(file.choose(), sep=",")


# Nome dos atributos
colnames(wine) <- c("Cvs", "Alcohol", "Malic acid", "Ash",
                    "Alcalinity of ash", "Magnesium", "Total phenols",
                    "Flavanoids", "Nonflavanoid phenols", "Proanthocyanins",
                    "Color intensity", "Hue", "OD280/OD315 of diluted wines",
                    "Proline")

# Classes
wine.classes <- factor(wine$Cvs)

# Calcula as componentes principais:
wine.pca <- prcomp(wine[,2:13], scale.=TRUE)

# Mostra resultados
summary(wine.pca)

# Carrega base:
data(swiss)

# Calcula e apresenta componentes principais:
swiss.pca <- prcomp(swiss[,1:6], scale.=TRUE)
swiss.pca

# Apresenta um resumo dos resultados:
summary(swiss.pca)

# Mostra as contribuicoes de cada componente principal:
screeplot(swiss.pca, type="barplot")


# Carrega base e calcula componentes principais:
pca <- prcomp(cars, scale.=TRUE)
summary(pca)

# Exemplo que ilustra o uso da matriz de covariancia e da matriz de correlacao
# na analise de componentes principais.
# Carregando a base
data(mtcars)
head(mtcars, 10)

# Calcula componentes principais com a matriz de covariancia
pc_cov <- prcomp(mtcars[,-11], scale=FALSE)

# Apresenta um resumo da analise
summary(pc_cov)

# Apresenta grafico com duas componentes principais
biplot(pc_cov)

# Calcula componentes principais com a matriz de correlacao
pc_cor <- prcomp(mtcars[,-11], scale=TRUE)

# Apresenta um resumo da analise
summary(pc_cor)

# Apresenta grafico com duas componentes principais
biplot(pc_cor)

# ICA no pacote R
# Carrega os pacotes necessarios e o conjunto de dados
library(fastICA)
data(iris)

# Aplica a tecnica ICA
iris.ica <- fastICA(iris[1:4], 2, maxit = 200, tol = 0.0001)

# Visualiza o resultado
plot(iris.ica$X, col=iris$Species, xlab="dimensão 1",
    ylab="dimensão 2", pch=16)

# IsoMap no pacote R
# Carrega os pacotes necessarios e o conjunto de dados
library(vegan)
data(iris)

# Aplica a tecnica IsoMap
iris.isomap <- isomap(vegdist(iris[1:4]), ndim=2, epsilon=0.3)

# Visualiza o resultado
plot(iris.isomap$points, col=iris$Species,
     xlab="dimensão 1", ylab="dimensão 2", pch=16)

# LLE no pacote R
# Carrega os pacotes necessarios e o conjunto de dados
source('auxiliar_aula3.R')
library(MASS)
data(iris)

# Aplica a tecnica LLE
iris.lle <- lle(iris[,1:4], m=2, k=2, reg=2,
                ss=FALSE, id=TRUE, v=0.9)

# Visualiza o resultado
plot(iris.lle$Y, col=iris$Species,
     xlab="dimensão 1", ylab="dimensão 2", pch=16)


# Sammon Mapping no pacote R
# Carrega os pacotes necessarios e o conjunto de dados
library(MASS)
data(iris)

# Remove dados duplicados (necessario para aplicar a tecnica)
iris_unicos <- unique(iris)

# Aplica a tecnica Sammon Mapping
iris.sammon <- sammon(dist(iris_unicos[,1:4]), k=2)

# Visualiza o resultado
plot(iris.sammon$points, col=iris_unicos$Species,
     xlab="dimensão 1", ylab="dimensão 2", pch=16)

# MDS no pacote R
# Carrega os pacotes necessarios e o conjunto de dados
library(MASS)
data(iris)

# Remove dados duplicados (necessario para aplicar a tecnica)
iris_unicos <- unique(iris)

# Aplica a tecnica MDS
iris.mds <- isoMDS(dist(iris_unicos[,1:4]), k=2)

# Visualiza o resultado
plot(iris.mds$points, col=iris_unicos$Species,
     xlab="dimensão 1", ylab="dimensão 2", pch=16)

# t-SNE no pacote R
# Carrega os pacotes necessarios e o conjunto de dados:
library(Rtsne)
library(datasets)
data(iris)

# Remove dados duplicados(necessario para executar t-SNE):
iris_unicos <- unique(iris)

# Executa t-SNE:
set.seed(42) # semente fixa para reprodutibilidade
tsne <- Rtsne(as.matrix(iris_unicos[,1:4]), perplexity = 30, dims=3)

# Visualiza o resultado:
plot(tsne$Y, col=iris_unicos$Species, xlab="dimensao 1", ylab="dimensao 2",
     pch=16)

# Carrega os pacotes necessarios e o conjunto de dados:
library(Rtsne)

# # Leitura arquivo da base wine versão online:
# train <- read.csv("https://www.ic.unicamp.br/~helio/MNIST/train.csv")

# Leitura arquivo da base MNIST com arquivo offline
# baixado do Moodle:
train <- read.csv(file.choose())


train$label <- as.factor(train$label)

# Executa t-SNE(nao ha dados duplicados):
set.seed(1) # para reprodutibilidade
tsne <- Rtsne(train[,-1], dims = 2, perplexity = 30, verbose = TRUE,
              max_iter = 500)

# Visualiza o resultado:
colors <- rainbow(length(unique(train$label)))
names(colors) <- unique(train$label)

plot(tsne$Y, t='n', main="", xlab="dimensao 1", ylab="dimensao 2")

text(tsne$Y, labels=train$label, col=colors[train$label], cex=0.5)


# UMAP no Pacote R
# Carrega os pacotes necessarios e o conjunto de dados:
library(umap)
data(iris)

# Executa UMAP:
set.seed(42) # semente fixa para reprodutibilidade
iris.umap <- umap(as.matrix(iris[,1:4]))

# Visualiza o resultado
plot(iris.umap$layout , col=iris$Species , xlab="dimensao 1",
     ylab="dimensao2", pch =16)

# Carrega os pacotes necessarios e o conjunto de dados:
library(umap)
# # Leitura arquivo da base wine versão online:
# train <- read.csv("https://www.ic.unicamp.br/~helio/MNIST/train.csv")

# Leitura arquivo da base MNIST com arquivo offline
# baixado do Moodle:
train <- read.csv(file.choose())

train$label <- as.factor(train$label)

# Executa UMAP:
set.seed(1) # para reprodutibilidade
mnist_umap <- umap(as.matrix(train[,-1]), n_neighbors=15,
                   min_dist=0.001, verbose=TRUE)

# Visualiza o resultado:
colors <- rainbow(length(unique(train$label)))
names(colors) <- unique(train$label)

plot(mnist_umap$layout , t='n', main="",
     xlab="dimensao 1", ylab="dimensao 2")
text(mnist_umap$layout, labels=train$label, col=colors[train$label], cex=0.5)


# Compara com a Analise de Componentes Principais:
pca <- princomp(train[,-1])$scores[,1:2]

plot(pca, t='n', main="", xlab="componente 1", ylab="componente 2")
text(pca, labels=train$label, col=colors[train$label])

##### Apêndice

# Considerando 2 atributos
# calcula vetor da media
apply(iris[,1:2], 2, mean)

# calcula matriz de covariancia
cov(iris[,1:2])

# Considerando 4 atributos
# calcula vetor da media
apply(iris[,1:4], 2, mean)
# calcula matriz de covariancia
cov(iris[,1:4])


# Autovalores e Autovetores
eig <- eigen(cov(iris[,1:4]))

# Escolha das Componentes Principais
eig <- eigen(cov(iris[,1:4]))
eig$vectors[,1:2]

# Construção da nova base após a redução
# prepara a matriz de dados
X <- iris[,1:4]

# centraliza os dados
X.scaled <- scale(X, center=TRUE)

# calcula matriz de covariancia
X.cov <- cov(X.scaled)
X.cov

# calcula autovalores/autovetores da matriz de correlacao
X.eig <- eigen(X.cov)
X.eig

# calcula a nova matriz de dados
X.new <- t(X.eig$vectors) %*% t(X.scaled)

# calcula a matriz transposta e renomeia colunas
X.new <- t(X.new)
colnames(X.new) <- c("PC1", "PC2", "PC3", "PC4")
head(X.new)


# Exemplo de redução de dimensionalidade
# cria matriz de entrada
X1 <- c(2.5, 0.5, 2.2, 1.9, 3.1, 2.3, 2.0, 1.0, 1.5, 1.1)
X2 <- c(2.4, 0.7, 2.9, 2.2, 3.0, 2.7, 1.6, 1.1, 1.6, 0.9)

X <- data.frame(X1, X2)
# centraliza dados
X.scaled <- scale(X, center=TRUE, scale=FALSE)
head(X.scaled)

# Grafico dos dados originais:
# visualiza caracteristicas
plot(X, xlab="X1", ylab="X2", col="blue", cex=1.0, pch=19)


# calcula matriz de covariancia
X.cov <- cov(X.scaled)

# calcula autovalores e autovetores
X.eig <- eigen(X.cov)
X.eig$values
X.eig$vectors

# Projecao dos dados com duas componentes principais:
# projeta dados no novo espaco
X.new <- t(X.eig$vectors) %*% t(X.scaled)
X.new <- t(X.new)

# Transformacao dos dados para o espaco original:
# transforma dados para o espaco original
X.back <- X.new %*% t(X.eig$vectors) + colMeans(X)


# Projecao dos dados com uma componente principal:
# projeta dados no novo espaco
X.new <- t(X.eig$vectors[,1]) %*% t(X.scaled)
X.new <- t(X.new)

# Transformacao dos dados para o espaco original:
# projeta dados no novo espaco
X.back <- X.new %*% X.eig$vectors[,1] + colMeans(X)
X.back


pizza <- read.csv("https://www.ic.unicamp.br/~helio/datasets/pizza.csv")
summary(pizza)

table(pizza$brand)

sapply(pizza, function(x) sum(is.na(x)))

pizza <- subset(pizza, select=-c(id))

pizza_cor <- cor(pizza[,2:length(pizza)])

library(corrplot)
corrplot(pizza_cor, method="square", type="upper",
         tl.col="black", tl.srt=45,
         addCoef.col="black", number.cex=0.7)

pca <- prcomp(pizza[,-1], center=TRUE, scale=TRUE)
summary(pca)

library(factoextra)
fviz_screeplot(pca, addlabels=TRUE, title="", x="Dimensções",
               y="Percentagem de variância explicada")

fviz_pca_var(pca, col.var="contrib",
             repel=TRUE, axes=c(1, 3)) +
   labs(title="", x="PC1", y="PC2")

fviz_pca_var(pca, col.var="contrib",
             repel=TRUE, axes=c(2, 3)) +
   labs(title="", x="PC2", y="PC3")

fviz_contrib(pca, choice="var", axes=1:2) +
   labs(title="", y="Contribuições (%)")


library(Rtsne)

set.seed(42) # semente fixa para reprodutibilidade
tsne <- Rtsne(as.matrix(pizza[,-1]), dims=3,
              check_duplicates=FALSE, perplexity=30)


colors <- rainbow(length(unique(pizza$brand)))
names(colors) <- unique(pizza$brand)

plot(tsne$Y, t='n',
     main="", xlab="dimensão 1", ylab="dimensão 2")

text(tsne$Y, labels=pizza$brand, col=colors[pizza$brand])


library(scatterplot3d)

scatterplot3d(x=tsne$Y[,1], y=tsne$Y[,2], z=tsne$Y[,3],
              color=colors[pizza$brand], main="",
              xlab="dimensão 1",
              ylab="dimensão 2",
              zlab="dimensão 3")

legend("topright", legend=c("A", "B", "C", "D", "E",
                            "F", "G", "H", "I", "J"),
       col=colors, xpd=TRUE, inset = -0.1,
       horiz=FALSE, ncol=2, pch=16)
