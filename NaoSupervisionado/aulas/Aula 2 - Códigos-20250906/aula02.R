##################################################################
# Mineracao de Dados Complexos -- MDC 
# Aprendizado de Máquina Não Supervisionado
# Prof. Helio Pedrini
# Codigos da Aula 2 - Regras de Associacao
##################################################################


# O pacote arules deve ser inicialmente instalado e carregado:
# install.packages("arulesCBA")
library(arulesCBA)

# Carrega base de dados
data(iris)
# Discretiza atributos continuos em intervalos usando modelo 
# de Fayyad e Irani:
set.seed(42) # para reprodutibilidade
dados <- discretizeDF.supervised(Species ~ ., data=iris)

# Apresenta dados
head(dados)
# Cria conjunto de treinamento com 90% dos dados
treino <- sample(1:nrow(dados), 
                 size = as.integer(nrow(dados)*0.9))

# Define valor de suporte
5/length(treino)

# Constroi um modelo CBA
modelo_cba <- CBA(Species ~ ., 
                  data = dados[treino,],
                  supp = 0.03, conf = 0.5)
modelo_cba

# Inspeciona regras
inspect(modelo_cba$rules)

# Realiza predicoes no conjunto de teste
predicao <- predict(modelo_cba, dados[-treino,])

# Apresenta os resultados
library(caret)
confusionMatrix(predicao, dados[-treino,]$Species)


# O pacote arules deve ser inicialmente instalado e carregado:
# install.packages("arules")
library(arules)

# Carrega os arquivos, format="basket" indica que o 
# arquivo esta no formato cesta de compras
transacoes <- read.transactions(file.choose(), format="basket", 
                                sep=",")

# As transacoes podem ser visualizadas usando a funcao inspect:
inspect(transacoes)

# As transacoes podem ser visualizadas com a funcao image:
image(transacoes)

# As regras de associacao podem ser produzidas com o uso da funcao
# apriori:
regras <- apriori(transacoes, parameter=list(supp=0.5, conf=0.5))

# Digitando o nome da variavel regras, tem-se o numero de regras
# criadas:
regras 

# O comando inspect permite a visualizacao das regras criadas:
inspect(regras)

# O pacote arulesViz possui varias funcionalidades de visualizacao 
# de regras
# install.packages("arulesViz")
library(arulesViz)

# Observacao: talvez seja necessario instalar o pacote git2r:
# install.packages("git2r")

# Para a visualizacao das regras, pode-se empregar a funcao plot:
plot(regras, method="graph") 

# Usando outro conjunto de dados (groceries.csv):
# carrega as bibliotecas
library(arules)
library(arulesViz)
library(datasets)

# carrega a base de dados
data(Groceries)
transacoes <- Groceries

# Resumo do conjunto de dados:
summary(transacoes)

# Histograma para os 20 primeiros itens:
itemFrequencyPlot(Groceries, topN = 20, type = "absolute")

# Histograma levando-se em conta a frequência do padrao na regra:
itemFrequencyPlot(transacoes, support = 0.1, cex.names = 0.8)

# Histograma levando-se em conta a frequência do padrao na regra.
# Com a reducao do valor do suporte, mais itens sao mostrados
# no grafico:
itemFrequencyPlot(transacoes, support=0.05, cex.names=0.8)

# Para a mineracao das regras, alguns limiares sao passados:
# suporte minimo igual a 0.001 e confianca minima igual a 0.8
regras <- apriori(Groceries, parameter = 
                    list(supp = 0.001, conf = 0.8))
regras

# Um resumo com as informacoes estatisticas das regras pode ser 
# obtido:
summary(regras)

# As cinco primeiras regras sao mostradas:
inspect(regras[1:5])

# As regras podem ser ordenadas de acordo com a relevância das
# medidas, por exemplo, a confianca:
regras <- sort(regras, by="confidence", decreasing=TRUE)

# As cinco primeiras regras em ordem de confianca sao:
inspect(regras[1:5])

# O que os clientes provavelmente comprariam antes de comprar
# leite integral?
regras <- apriori(data = Groceries, 
                  parameter=list(supp=0.001, conf=0.08), 
                  appearance=list(default="lhs", rhs="whole milk"), 
                  control=list(verbose=F))
regras <- sort(regras, decreasing=TRUE, by="confidence")
inspect(regras[1:5])

# O que os clientes provavelmente comprariam se comprassem leite
# integral?
regras <- apriori(data=Groceries, 
                  parameter=list(supp=0.001, conf=0.15,minlen=2), 
                  appearance=list(default="rhs", lhs="whole milk"),
                  control=list(verbose=F))
regras <- sort(regras, decreasing=TRUE, by="confidence")
inspect(regras[1:5])

# Para a visualizacao das regras:
plot(regras, method="graph") 

# Mineracao de regras com outros valores de suporte e confianca:
regras <- apriori(Groceries, parameter=list(supp=0.009,
                  conf=0.25, minlen=2))
regras

# Analise das primeiras cinco regras ordenadas por valor de lift:
inspect(head(sort(regras, by="lift"), 5))

# Analise das regras com alto suporte e alta confianca:
inspect(sort(sort(regras, by="support"), by="confidence")[1:5])

# Pode-se exibir as regras em que os eixos representam o suporte e 
# a confianca, enquanto as cores representam valores de lift:
plot(regras, measure=c("support", "confidence"), shading="lift", 
     jitter=0)


# Pode-se selecionar um produto e verificar a quais produtos ele
# esta associado ou se associa:
meat.regras <- sort(subset(regras, 
                           subset=lhs %in% "beef" | 
                             lhs %in% "sausage" | 
                             lhs %in% "chicken"),
                    by="confidence")
summary(meat.regras)

# Visualizacao das regras criadas:
inspect(meat.regras)

# A visualizacao das regras para o produto meat pode ser 
# realizada como:
plot(meat.regras, method="graph",  shading="lift")

# Outros produtos podem ser inspecionados individualmente:
milk.regras <- sort(subset(regras, subset = rhs %in% "whole milk"), 
                    by = "confidence")
summary(milk.regras)
inspect(milk.regras)

coke.regras <- sort(subset(regras, subset = rhs %in% "soda"), 
                    by = "confidence")
summary(coke.regras)
inspect(coke.regras)

yogurt.regras <- sort(subset(regras, subset = lhs %in%  "yogurt"),  
                      by = "confidence")
summary(yogurt.regras)
inspect(yogurt.regras)

# As regras para o produto meat podem ser exibidas na forma 
# matricial:
library(colorspace)
plot(meat.regras, method="matrix",
     measure=c("support", "confidence"), 
     control=list(col=sequential_hcl(3)))

# Uma forma conveniente e o grafico de matriz agrupada:
plot(meat.regras, method="grouped", measure="support", 
     control=list(col=sequential_hcl(3)))

# Outro grafico interessante para mostrar quantos produtos 
# de cada tipo estao disponiveis e chamado de mapa de arvores
# install.packages("treemap")
# install.packages("magrittr")
# install.packages("dplyr")
library(treemap)
library(magrittr)
library(dplyr)


# Mapa de arvores com nivel de segmentacao para apresentar quais
# produtos estao nas listas de compras
occur1 <- transacoes@itemInfo %>% 
          group_by(level1) %>% 
          summarize(n = n())
treemap(occur1,index = c("level1"), vSize = "n", 
        title = "", palette = "Dark2", border.col = "#FFFFFF")

# Mapa de arvores com nivel maior de segmentacao para apresentar quais
# produtos estao nas listas de compras
occur2 <- transacoes@itemInfo %>% 
          group_by(level1, level2) %>% 
          summarize(n = n())
treemap(occur2, index = c("level1",  "level2"), vSize = "n", 
        title = "", palette = "Dark2", border.col = "#FFFFFF")

# Mapa de arvores com nivel ainda maior de segmentacao para 
# apresentar quais produtos estao nas listas de compras
occur3 <- transacoes@itemInfo %>% 
          group_by(level1, level2, labels) %>% 
          summarize(n = n())
treemap(occur3, index = c("level1", "labels"), vSize = "n", 
        title = "", palette = "Dark2", border.col = "#FFFFFF")


# Um teste interessante e identificar quais itens sao menos 
# provaveis de serem comprados juntos. Isto pode ser feito 
# verificando se o valor de lift e menor do que 1:
inspect(tail(sort(regras, by="lift")))

# Calculo de varias medidas de interesse:
regras <- apriori(Groceries, 
                  parameter = list(supp = 0.01, 
                                   conf = 0.25, 
                                   minlen = 2))


medidas <- interestMeasure(regras, 
                           c("confidence", "conviction",
                             "coverage", "support", 
                             "oddsRatio", "leverage"), Groceries)

head(medidas)

# Carrega pacotes e base de dados:
library(arules)
library(arulesViz)
data("AdultUCI")

# Apresenta informacoes basicas sobre o conjunto de dados:
dim(AdultUCI)

# Apresenta seis primeiras linhas da base:
head(AdultUCI)

# Remove dois atributos:
AdultUCI[["fnlwgt"]] <- NULL
AdultUCI[["education-num"]] <- NULL

# Converte quatro atributos contínuos para ordinais:
AdultUCI[["age"]] <-  ordered(cut(AdultUCI[["age"]], 
                              c(15,25,45,65,100)), labels = c("Young", 
                               "Middle-aged", "Senior","Old"))
AdultUCI[["hours-per-week"]] <- 
              ordered(cut(AdultUCI[["hours-per-week"]],
                          c(0,25,40,60,168)), 
                      labels = c("Part-time",
                                 "Full-time",
                                 "Over-time",
                                 "Workaholic"))
AdultUCI[["capital-gain"]] <- 
    ordered(cut(AdultUCI[["capital-gain"]],
            c(-Inf, 0, 
              median(AdultUCI[["capital-gain"]]
                     [AdultUCI[["capital-gain"]] > 0]),
              Inf)),
            labels = c("None", "Low", "High"))
AdultUCI[["capital-loss"]] <- 
  ordered(cut(AdultUCI[["capital-loss"]],
              c(-Inf, 0,
                median(AdultUCI[["capital-loss"]]
                       [AdultUCI[["capital-loss"]] > 0]),
                Inf)),
              labels = c("None", "Low", "High"))

# Converte dados para formato de transacoes:
Adult <- as(AdultUCI, "transactions")

# Apresenta transacoes:
Adult
summary(Adult)

# Mostra as frequências de itens no conjunto de dados com suporte
  # maior do 10%:
itemFrequencyPlot(Adult, support=0.10)

# Ha quantas regras com suporte mı́nimo de 1% e confianca de 60%?
# O comando a seguir exibe o warning sobre a quantidade maxima 
# de itens em uma regra
regras <- apriori(Adult, parameter=list(
                support=0.01, confidence=0.6, minlen=2))

# Resumo das regras:
summary(regras)

# Busca todas as regras relacionadas a variavel 
# income (high income e low income):
regras_high_income <- apriori(Adult, parameter=list(
      supp=0.01, conf=0.6, minlen=2), appearance = list(
      default="rhs", lhs="income=large"))

# Resumo das regras relacionadas a variavel income:
summary(regras_high_income)

# Regras para low income:
regras_low_income <- apriori(Adult, parameter=list(
      supp=0.01, conf=0.6, minlen=2), appearance = list(
      default="rhs", lhs="income=small"))

# Resumo das regras relacionadas a variavel income:
summary(regras_low_income)

# Se houver regras redundantes, poda-las:
summary(is.redundant(regras_high_income))
regras_high_income[!is.redundant(regras_high_income)]

# Se houver regras redundantes, poda-las:
summary(is.redundant(regras_low_income))
regras_low_income[!is.redundant(regras_low_income)]

# Inspeciona os três conjuntos de itens com maior lift que 
# indicam low income:
inspect(head(sort(regras_low_income, by ="lift"), 3))

# Inspeciona os três conjuntos de itens com maior lift que 
# indicam high income:
inspect(head(sort(regras_high_income, by ="lift"), 3))

# Mostra as regras para low income suporte e confianca:
plot(regras_low_income, method="graph")

# Mostra as regras para high income suporte e confianca:
plot(regras_high_income, method="graph")
