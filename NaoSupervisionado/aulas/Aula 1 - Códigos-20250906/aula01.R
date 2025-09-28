##################################################################
# Mineracao de Dados Complexos -- MDC 
# Aprendizado de Máquina Não Supervisionado
# Prof. Helio Pedrini
# Codigos da Aula 1 - Conceitos, Tecnicas e Aplicacoes
##################################################################

### Matriz de Confusao em R
# Carregando bibliotecas
library(caret)
library(e1071)
# Exemplo: matriz de confusao com 2 classes
verdade <- factor(c(1, 1, 0, 1, 0, 0, 1, 0, 0, 0))
predito <- factor(c(1, 0, 0, 1, 0, 0, 1, 1, 1, 0))

resultado <- confusionMatrix(predito, verdade)
print(resultado)

# Exemplo: matriz de confusao com 2 classes
classes <- c("normal", "abnormal")
verdade <- factor(rep(classes, times = c(86, 258)),
                  levels = rev(classes))
predito <- factor(c(rep(classes, times = c(54, 32)),
                    rep(classes, times = c(27, 231))),
                  levels = rev(classes))

tabela <- table(predito, verdade)
confusionMatrix(tabela)
confusionMatrix(predito, verdade)

# Exemplo: matriz com 3 classes
verdade <- factor(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2,
                    2, 2, 2, 2, 2, 2, 2, 2, 2, 2))
predito <- factor(c(0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0,
                    0, 0, 0, 1, 1, 1, 1, 2, 1, 1, 2,
                    2, 2, 2, 2, 2, 2, 2, 2, 2, 2))

tabela <- table(predito, verdade)
resultado <- confusionMatrix(tabela)
print(resultado)

# Exemplo: matriz de confusao com 3 classes
confusionMatrix(iris$Species, sample(iris$Species))

newPrior <- c(.05, .8, .15)
names(newPrior)<- levels(iris$Species)

confusionMatrix(iris$Species, sample(iris$Species))

# Exemplo 4: matriz de confusao com 3 classes (base iris)

## Lembrete: as 4 primeiras colunas da base iris contem as
# medidas de altura e largura de sepala e petala
TrainData <- iris[,1:4]

## Lembrete: a 5a coluna da base iris indica
# a especie da amostra
TrainClasses <- iris[,5]

knnFit <- train(TrainData, TrainClasses, method = "knn",
                preProcess = c("center", "scale"),
                tuneLength = 10,
                trControl = trainControl(method = "cv"))

confusionMatrix(knnFit)
confusionMatrix(knnFit, "average")
confusionMatrix(knnFit, "none")


#### Analise exploratoria em R
## Base faithful: colecao de observacoes do geiser
X <- faithful

# Apresenta as primeiras seis linhas:
head(X)

# Informacoes basicas do dataframe:
# numero de linhas:
nrow(X)
# numero de colunas:
ncol(X)
# dimensoes:
dim(X)
# nomes das colunas:
names(X)

# Media: valor central do conjunto de dados
duracao <- X$eruptions
mean(duracao)

# Variancia: mede quao afastados da media estao os
# valores do conjunto.
duracao <- X$eruptions
var(duracao)

# Desvio Padrao: medida de dispersao que corresponde a
# raiz quadrada da variancia.
duracao <- X$eruptions
sd(duracao)

# Covariancia: mede se duas variaveis x e y estao linearmente
# relacionadas.
duracao <- X$eruptions
espera <- X$waiting
cov(duracao, espera)

# Coeficiente de Correlacao: medida normalizada de quao
# relacionadas linearmente duas variaveis estao.
duracao <- X$eruptions
espera <- X$waiting
cor(duracao, espera)

# Momento Central: o k-esimo momento central de um conjunto.
library(e1071)
duracao <- X$eruptions
moment(duracao, order = 3, center = TRUE)

# Assimetria ou Obliquidade: mede a falta de simetria de uma
# determinada distribuicao de frequencia.
library(e1071)
duracao <- X$eruptions
skewness(duracao)

# Curtose: descreve a forma da distribuicao dos dados, que
# caracteriza o achatamento da funcao de distribuicao
# de probabilidade.
library(e1071)
duracao <- X$eruptions
kurtosis(duracao)

# Minimo: menor valor observado.
duracao <- X$eruptions
min(duracao)

espera <- X$waiting
min(espera)

# Maximo: maior valor observado.
duracao <- X$eruptions
max(duracao)

espera <- X$waiting
max(espera)

# Amplitude: medida simples de dispersao dos dados expressa
# pela diferenca entre o valor maximo e o valor
# minimo dos dados.
duracao = X$eruptions
max(duracao) - min(duracao)

espera <- X$waiting
max(espera) - min(espera)


# Mediana: valor que separa a metade maior e menor
# do conjunto ordenado de dados.
duracao <- X$eruptions
median(duracao)

# Moda: valor mais frequente no conjunto de observacoes.
moda <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match (v, uniqv)))]
}

moda(X$eruptions)
moda(X$waiting)

# Quartis: valores que dividem o conjunto ordenado de dados
# em quatro  partes iguais, tal que cada parte
# representa 1/4 das observacoes
duracao <- X$eruptions
quantile(duracao)

# Percentis: valores que dividem o conjunto ordenado de
# dados em 100  partes, cada um com uma percentagem de dados
# aproximadamente igual.
duracao <- X$eruptions
quantile(duracao, c (.32, .57, .98))

# Intervalo Interquartil: grau de espalhamento (dispersao)
# dos dados em torno da medida de centralidade.
duracao <- X$eruptions
IQR(duracao)


# Box Plot: representacao grafica baseada nos quartis da
# variavel observada.
duracao <- X$eruptions
boxplot(duracao, horizontal = TRUE)

# Scatter Plot: grafico de dispersao com pontos geometricos.
duracao <- X$eruptions
espera <- X$waiting
plot(duracao, espera,
     xlab = "Duracao da Erupcao",
     ylab = "Tempo Esperado")


## Base iris
# Apresenta um resumo da base:
summary(iris)

# Apresenta proporcao de amostras para cada classe:
tabela <- table(iris$Species)
rotulos <- paste(names(tabela), "\n", tabela, sep="")
pie(tabela, labels=rotulos)

# Apresenta distribuicao de atributos entre as classes:
library(reshape2)
iris_melted <- melt(iris, id="Species")

library(ggplot2)
bar <- ggplot(data=iris_melted, aes(x=Species, y=value,
                                    fill=variable))
bar + geom_bar(stat="identity", position="dodge") + 
  scale_fill_manual(values=c("orange", "blue",
                             "darkgreen", "purple"),
                    name="Iris\nMeasurements",
                    breaks=c("Sepal.Length", "Sepal.Width",
                             "Petal.Length", "Petal.Width"),
                    labels=c("Sepal Length", "Sepal Width",
                             "Petal Length", "Petal Width"))

# Histograma circular de atributos por classe:
library(ggplot2)
ggplot(iris, aes(x=Sepal.Length, fill=Species)) +
  geom_bar(width=0.1, color="white") + coord_polar()

# Instala e carrega o pacote dplyr para manipular dados:
## install.packages("dplyr")
library(dplyr)

# Filtra dados para a especie virginica:
names(iris) <- tolower(names(iris))
virginica <- filter(iris, species == "virginica")

# Apresenta as primeiras seis linhas:
head(virginica)

# Filtra com multiplas condicoes:
sepalLength6 <- filter(iris, species == "virginica",
                       sepal.length > 6)

# Apresenta as ultimas seis linhas:
tail(sepalLength6)

# Seleciona colunas especificas:
selected1 <- dplyr::select(iris, sepal.length,
                           sepal.width, petal.length)

# Seleciona intervalo de colunas:
selected2 <- dplyr::select(iris, 
                           sepal.length : petal.length)

# Apresenta as primeiras tres linhas:
head(selected1, 3)

# Verifica resultados selecionados:
identical(selected1, selected2)

# Outra forma de selecionar amostras (5 primeiras linhas):
subset(iris, species == "setosa")[1:5,]
subset(iris, species == "versicolor")[1:5,]
subset(iris, species == "virginica")[1:5,]

# Visualiza os dados:
plot (iris)

# Visualiza os dados conforme ocorrencia de certo atributo:
plot(iris$sepal.width, iris$sepal.length,
     xlab = "sepal width", ylab = "sepal length")

# Visualiza distribuicao
hist(iris$sepal.width,main = "", xlab = "sepal width", 
     ylab = "frequency")

hist(iris$petal.length,main="", xlab="petal length", 
     ylab="frequency")

# Visualiza estatísticas
par(mar = c(7, 5, 1, 1)) #aumenta o espaco para os rotulos
boxplot(iris,las=2)

irisVer <- subset(iris,species=="versicolor")
irisSet <- subset(iris,species=="setosa")
irisVir <- subset(iris,species=="virginica")

par (mfrow = c(1, 3),mar = c(6, 3, 2, 1))
boxplot (irisVer[,1:4], main = "versicolor", ylim=c(0, 8), 
         las = 2)
boxplot(irisSet[,1:4],main = "setosa", ylim = c(0, 8), 
        las = 2)
boxplot(irisVir[,1:4], main = "virginica", ylim = c(0, 8), 
        las = 2)

# Apresenta distribuicao de dados para atributo e classe
par(mfrow = c(1, 3))
hist(irisVer$petal.length, breaks = seq(0, 8, l = 17), 
     main = "", xlab = "petal length", ylab = "frequency", 
     xlim = c(0, 8), ylim = c(0, 40))
hist(irisSet$petal.length, breaks = seq(0, 8, l = 17), 
     main = "", xlab = "petal length",ylab = "frequency", 
     xlim =c(0, 8), ylim = c(0, 40))
hist(irisVir$petal.length,breaks = seq(0, 8, l = 17), 
     main = "", xlab = "petal length",ylab = "frequency", 
     xlim = c(0, 8), ylim = c(0, 40))

# Correlacao entre variaveis
corr <- cor(iris[,1:4])
round(corr, 3)

# Visualiza correlacao entre classes
pairs (iris[,1:4], col = iris[,5], oma = c(4, 4, 6, 12))
par(xpd = TRUE)
legend(0.85, 0.6, as.vector(unique(iris$species)), 
       fill = c(1, 2, 3))

# Visualiza por coordenadas paralelas
library(MASS)
par(mfrow = c(1, 1))

parcoord(iris[,1:4], col = iris[,5],var.label = TRUE,
         oma = c(4, 4, 6, 12))
par(xpd = TRUE)

legend(2.2, 1.1, as.vector(unique(iris$species)),fill = c(1, 2, 3))

## Exemplo base de dados cars
# Lista classes
names(cars)

# Apresenta resumo da base 
summary(cars)

# Primeiros 6 registros
head(cars)

# Ultimos 6 registros
tail(cars)

# Verifica existencia de dados faltantes
sapply(cars, function(x) sum(is.na(x)))

# Visualiza relacao entre dados
plot(cars$speed, cars$dist, main = "", col = "red", 
     xlab = "speed", ylab = "dist")

out = lm(dist ~ speed, data=cars)
plot(cars)
abline(out, col="red")

# Verifica existencia de outliers
par(mfrow = c(1, 2))

boxplot(cars$speed, main = "speed",
        sub = paste("outlier rows : ",
                    boxplot.stats(cars$speed)$out))
boxplot(cars$dist, main = "dist", 
        sub = paste("outlier rows :", 
                    boxplot.stats(cars$dist)$out))

# Histograma dos atributos
par(mfrow = c(1, 2)) 

hist(x = cars$speed, col = "red", main = "", 
     xlab = "Velocidade (mph)", ylab = "Frequência")
hist(x = cars$dist, col = "lightblue", main = "", 
     xlab = "Distância (ft)", ylab = "Frequência")

#recomocao de dados discrepantes
outliers <- boxplot(cars$dist, plot = FALSE)$out
print(outliers)
cars <- cars[-which(cars$dist %in% outliers),]
par(mfrow = c(1, 1))
boxplot(cars$dist, plot = TRUE)

## Exemplo dados duplicados e faltantes
# Carrega pacote
library(dplyr)

# Cria vetor
x <- c(1, 1, 4, 5, 4, 6)

# Encontra posicao dos elementos duplicados
duplicated(x)

# Obtem elementos duplicados
x[duplicated(x)]

# Remove elementos duplicados
x[!duplicated(x)]

# Obtem elementos unicos
unique(x)

# Cria vetor
y <- c(1, 4, 16, -5, 25, -6)

# Operacao que cria dados inválidos
z <- sqrt(y)
z

# Encontra posicoes faltantes
is.na(z)
which(is.na(z))

# Realiza operacao com dados faltantes (incorreto)
mean(z)

# Realiza operacao com dados faltantes (correto)
mean(z,na.rm = TRUE)

#preenche dados faltantes com a media
z[is.na(z)] <- mean(z,na.rm =TRUE)
z

# Tabela com dados faltantes
df <- data.frame(col1 = c(1:3, NA),
                 col2 = c( "parte", NA, "de", "texto"), 
                 col3 = c(TRUE, FALSE, TRUE, TRUE), 
                 col4 = c(2.5, 4.2, 3.2, NA), 
                 stringsAsFactors = FALSE)

# Encontra linhas completas
complete.cases(df)

# Subconjunto de dados completos
df[complete.cases(df),]

# Subconjunto de dados faltantes
df[!complete.cases(df),]

# Omite linhas faltantes
na.omit(df)


## Exemplo conversao de dados e normalizacao
# Estrategia one hot encoding
df <- data.frame(Marca = c("Volkswagen", "Toyota", 
                           "Honda", "Honda"), 
                 Preco = c(30000, 35000, 37000, 41000),
                 stringsAsFactors = FALSE)

library(caret)
zz = dummyVars(~.,df)
df_novo <- predict(zz,df)
df_novo

# Outra forma de aplicar a conversao
library(mltools)
library(data.table)

df_novo <- one_hot(as.data.table(df))


# Normalizacao de variaveis numericas
# Cria tabela de dados
idade <- c(25, 35, 40, 50)
salario <- c(200000, 1200000, 160000, 2000000)

df <- data.frame("Idade" = idade, "Salario" = salario,
                 stringsAsFactors = FALSE)
df


# Normalizacao z-score
df_zscore <- as.data.frame(scale(df))
df_zscore

# Normalizacao min-max
normalize <- function(x){
  return((x - min(x)) / (max(x) - min(x)))
}

df_minmax <- as.data.frame(lapply(df,normalize))
df_minmax

# Tratamento de dados desbalanceados
source('auxiliar_aula1.R')

# le dados de entrada
# retorna para a base original, onde os nomes dos atributos
# comecam com letra maiuscula
data(iris) 
data <- iris[,c(1, 2, 5)]

# altera dados
data$Species <- factor(ifelse(data$Species == "setosa",
                              "rare", "common"))

# verifica distribuicao das classes dos dados gerados
table(data$Species)

# utilizando SMOTE
new_data<- SMOTE(Species~., data, k = 3, perc.under= 100, 
                 perc.over = 600)

# verifica distribuicao
table(new_data$Species)



# Teste de normalidade 
# Exemplo 1: distribuicao aproximadamente normal
# (resultados podem ser ligeiramente diferentes dos slides
# por causa da randomizacao)
normal <- rnorm(100, mean = 5, sd = 3)
hist(normal, probability = T, main = "", ylab = "frequencia",
     xlab="dados com distribuicao aproximadamente normal")

# Teste de Shapiro-Wilk
shapiro.test(normal)

# Teste de Cramer-von Mises
library(nortest)
cvm.test(normal)

# Teste de Kolmogorov-Smirnov
lillie.test(normal)


# Teste de Anderson-Darling 
ad.test(normal)

# Exemplo: distribuicao assimetrica
skewed <- runif(100, min=2, max=4)
hist(skewed, probability=T, main="", ylab="frequência", 
     xlab="dados com distribuicao assimétrica")

# Teste de Shapiro-Wilk
shapiro.test(skewed) 

# Teste de Cramer-von Mises
cvm.test(skewed)

# Teste de Kolmogorov-Smirnov
lillie.test(skewed)

# Teste de Anderson-Darling 
ad.test(skewed)

# Deteccao de Anomalias
# Carrega a base de dados
library(ggplot2)
data(mpg)
head(mpg)
# Visualiza estrutura da base
str(mpg)

# Estatisticas Descritivas
# Minimo e maximo
min(mpg$hwy)
max(mpg$hwy)

# Resumo das estatisticas
summary(mpg$hwy)

# Histograma
hist(mpg$hwy, xlab="hwy", breaks=sqrt(nrow(mpg)))

# Boxplot
boxplot(mpg$hwy, ylab="hwy")

# Regra de tres desvios-padrao
# Limites inferior e superior
lower_bound <- mean(mpg$hwy) - 3 * sd(mpg$hwy)
lower_bound

upper_bound <- mean(mpg$hwy) + 3 * sd(mpg$hwy)
upper_bound

# Todas as observacoes do limite inferior e acima do limite superior sao
# consideradas possiveis outliers:
out_ind <- which(mpg$hwy < lower_bound | mpg$hwy > upper_bound)
out_ind
# Linhas correspondentes aos outliers
mpg[out_ind,]

# Intervalo Interquartil
# Calculo deo intervalo
Q <- quantile(mpg$hwy, probs=c(.25, .75))
iqr <- IQR(mpg$hwy)
lower_bound <- Q[1] - 1.5 * iqr
lower_bound

upper_bound <- Q[2] + 1.5 * iqr
upper_bound

# Linhas correspondentes aos outliers
out_ind <- which(mpg$hwy < lower_bound | mpg$hwy > upper_bound)
out_ind

# Pode-se tambem usar o boxplot para extrair os valores dos possiveis outliers
# pelo criterio do IQR:
boxplot.stats(mpg$hwy)$out


# Linhas correspondentes aos outliers
out <- boxplot.stats(mpg$hwy)$out
out_ind <- which(mpg$hwy %in% c(out))
out_ind

mpg[out_ind,]


# Percentis
# observacoes que estao fora do intervalo de 2.5 a 97.5 percentis
# serao consideradas possiveis outliers:
lower_bound <- quantile(mpg$hwy, 0.025)
lower_bound
upper_bound <- quantile(mpg$hwy, 0.975)
upper_bound

# Todas as observacoes abaixo de 14 e acima de 35.175 serao consideradas
# possiveis outliers
out_ind <- which(mpg$hwy < lower_bound | mpg$hwy > upper_bound)
out_ind
mpg[out_ind,]

# Atribuindo 1 e 99 aos percentis
lower_bound <- quantile(mpg$hwy, 0.01)
upper_bound <- quantile(mpg$hwy, 0.99)
out_ind <- which(mpg$hwy < lower_bound | mpg$hwy > upper_bound)
mpg[out_ind,]

# Filto de Hampel
lower_bound <- median(mpg$hwy) - 3 * mad(mpg$hwy)
lower_bound
upper_bound <- median(mpg$hwy) + 3 * mad(mpg$hwy)
upper_bound

# Linhas correspondentes aos outliers (0 linhas neste exemplo)
out_ind <- which(mpg$hwy < lower_bound | mpg$hwy > upper_bound)
out_ind

# Teste de Grubbs
library(outliers)

# Para testar o maior valor
test <- grubbs.test(mpg$hwy)
test

# Para testar o menor valor
test <- grubbs.test(mpg$hwy, opposite=TRUE)
test


# Teste de Dixon
submpg <- mpg[1:20,]
test <- dixon.test(submpg$hwy)
test

# Para testar o maior valor
test <- dixon.test(submpg$hwy, opposite=TRUE)
test

# comparando o teste estatistico com o boxplot
out <- boxplot.stats(submpg$hwy)$out
boxplot(submpg$hwy, ylab="hwy")
mtext(paste("Outliers: ", paste(out, collapse=", ")))

remove_ind <- which.min(submpg$hwy)
subsubmpg <- submpg[-remove_ind,]

test <- dixon.test(subsubmpg$hwy)
test

# Teste de Rosner
library(EnvStats)
test <- rosnerTest(mpg$hwy, k=3)
test

# Local Outlier Factor (LOF)
source('auxiliar_aula1.R')
outlier.scores <- lofactor(mpg$hwy, k=3)
# Seleciona os primeiros 3 potenciais outliers
outliers <- order(outlier.scores, decreasing=T)[1:3]
print(outliers)
# Visualiza amostras que possuem fator de densidade superior a 1.5
print(which(outlier.scores > 1.5))


# Normalizacao de cadeias de caracteres
library(stringr)

# Operacoes de limpeza de dados
# remove espacos no inicio e no fim da palavra
str_trim("  exemplo de texto ")
# remove espacos apenas no inicio da palavra
str_trim("  exemplo de texto ", side="left")
# remove espacos apenas no fim da palavra
str_trim("  exemplo de texto ", side="right")

# Conversao de cadeias de caracteres
# converte para letras maiusculas
toupper("Exemplo de texto") 
# converte para letras minusculas
tolower("Exemplo de texto")

# manipulacao de datas
library(lubridate)
x <- "09/01/02"

# conversao entre formatos de datas
# ano-mes-dia
ymd(x)
# mes-dia-ano
mdy(x)
# ano-dia-mes
ydm(x)
# mes-ano-dia
myd(x)
# dia-mes-y
dmy(x)
# dia-ano-mes
dym(x)

# obtem horario corrente
now()

now("GMT")


## Mineracao de Textos
# Instalar pacotes
# install.packages("tm")
# install.packages("RColorBrewer")
# install.packages("wordcloud")
# install.packages("ggplot2")
# install.packages("dplyr")

# Carregar pacotes
library(tm)
library(RColorBrewer)
library(wordcloud)
library(ggplot2)
library(dplyr)

# # Ler texto de entrada pelo link
# arquivo <- "http://www.sthda.com/sthda/RDoc/example-files/martin-luther-king-i-have-a-dream-speech.txt"
# texto <- readLines(arquivo)

# Ler texto de entrada local (escolhendo o arquivo mlk.txt
# disponibilizado no moodle pelo sistema de arquivos)
texto <- readLines(file.choose())

# Converte texto para formato corpus:
corpus <- Corpus(VectorSource(texto))

# Limpa corpus
corpus <- tm_map(corpus, stripWhitespace)
corpus <- tm_map(corpus, tolower)
corpus <- tm_map(corpus, removeNumbers)
corpus <- tm_map(corpus, removePunctuation)
corpus <- tm_map(corpus, removeWords, stopwords("english"))

# Mostra nuvem de palavras (max.words: numero maximo de 
# palavras que serao apresentadas):
wordcloud (corpus, scale=c(5, 0.5), max.words=200, 
           random.order=FALSE, rot.per=0.35, 
           use.r.layout=FALSE, colors=brewer.pal(8, "Dark2"))

# Mostra nuvem de palavras (palavras com frequencia abaixo 
# de min.freq nao serao apresentadas):
wordcloud(corpus, scale=c(5, 0.5), min.freq=1, 
          max.words=200, random.order=FALSE, rot.per=0.35, 
          use.r.layout=FALSE, colors=brewer.pal(8, "Dark2"))

# Converte corpus em matriz de documento de termos, em que 
# cada linha e uma palavra e cada documento e uma coluna:
tdm <- TermDocumentMatrix(corpus)

# Calcula frequencia de palavras
contagem <- apply(tdm, 1, sum)
freq <- contagem[order(contagem, decreasing=TRUE)]
head(freq, 20)

# Mostra histograma de palavras
freq_hist <- data.frame(word=names(freq), freq=freq)
subset(freq_hist, freq > 5) %>% ggplot(aes(word, freq)) + 
  geom_bar(stat = "identity", fill = "darkred", 
           colour = "darkgreen") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "palavras", y = "frequencia")

# Ao verificar o conteudo da matriz de termos, pode-se 
# identificar a esparsidade das palavras:
tdm

# Pode-se remover do corpus as palavras menos frequentemente
# utilizadas matriz deve ter no maximo 85% de espaco vazio
tdms <- removeSparseTerms(tdm, 0.85)
tdms

# Calcula frequ?encia de palavras
contagem <- apply(tdms, 1, sum)
freq <- contagem[order(contagem, decreasing=TRUE)]

# Mostra histograma de palavras
freq_hist <- data.frame(word=names(freq), freq=freq)
subset(freq_hist, freq > 5) %>% ggplot(aes(word, freq)) + 
  geom_bar(stat = "identity", fill = "darkred", 
           colour = "darkgreen") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "palavras", y = "frequencia")

# Dado um determinado termo, pode-se identificar as palavras 
# que sao mais altamente correlacionadas com ele: limite da 
# correlacao igual a 0.80
findAssocs(tdm, c("freedom"), corlimit = 0.80)

# Carregar o pacote
library(textstem)

# A lematizacao reduz as palavras flexionadas a forma raiz
# para agrupar variacoes morfologicas
texto <- c("fly", "flies", "flying")
lemmatize_words(texto)

# Carregando os pacotes
library(syuzhet)
library(tidytext)
library(dplyr)

texto <- "I love programming in R, but sometimes it is difficult."

# obtem a pontuacao de sentimento do texto com o lexico Bing
sentimento_texto <- get_sentiment(texto, method="bing")

# exibe a pontuacaoo geral do texto
print(sentimento_texto)

# tokeniza o texto em palavras e converte para um data frame
palavras <- data.frame(word=get_tokens(texto), 
                       stringsAsFactors=FALSE)

# obtem o lexico Bing (palavras e suas classificacoes)
sentimentos <- get_sentiments("bing")

# une as palavras com o lexico
sentimento_palavra <- merge(palavras, sentimentos, 
                            by="word", all.x=FALSE)

# exibe as palavras que contribuiram para o sentimento do texto
print(sentimento_palavra)

# carregando os pacotes
library(tidytext)
library(dplyr)

texto <- data.frame(text=c("Várias pessoas gostam de estudar 
                            matemática e gostam de programação."))

# captura padroes entre palavras por meio de n-grams
texto %>%
  unnest_tokens(ngram, text, token="ngrams", n=2) %>%
  count(ngram, sort=TRUE)

# carregando os pacotes
library(tidytext)
library(dplyr)

# exemplo de corpus com diferentes discursos
discursos <- data.frame(
  autor = c("Autor A", "Autor A", "Autor B", "Autor B"),
  texto = c("A economia está crescendo rapidamente",
            "O governo precisa investir em educação e tecnologia",
            "A tecnologia está mudando o mundo",
            "A inovação científica é fundamental")
)

# tokeniza palavras
tokens <- discursos %>%
  unnest_tokens(word, texto) %>%
  count(autor, word, sort=TRUE) %>%
  bind_tf_idf(word, autor, n)

# exibe palavras mais representativas por autor
tokens %>%
  group_by(autor) %>%
  slice_max(tf_idf, n=5) %>%
  ungroup()

