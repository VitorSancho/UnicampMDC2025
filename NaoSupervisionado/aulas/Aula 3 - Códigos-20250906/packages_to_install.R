##################################################################
# Mineração de Dados Complexos -- MDC 2025
# Aprendizado de Maquina Nao Supervisionado
#
# Abra este arquivo com o Rstudio e execute cada linha
# separadamente. Caso encontre algum erro entre em contato com os
# monitores. Lembre-se de indicar aos monitores o seu sistemaa
# operacional e a versão do R instalada.
##################################################################


##################################################################
# Pacotes para a Aula 3
##################################################################

# Instalando o pacote
install.packages("readr")
library(readr)
# Mensagem esperada:
# (Sem mensagem)

# Instalando o pacote
install.packages("nnls")
# Carregando o pacote
library(nnls)
# Mensagem esperada:
# (Sem mensagem)

# Instalando o pacote
install.packages("rsvd")
# Carregando o pacote
library(rsvd)
# Mensagem esperada:
# (Sem mensagem)

# Instalando o pacote
install.packages("mlbench")
# Carregando o pacote
library(mlbench)
# Mensagem esperada:
# (Sem mensagem)

# Instalando o pacote
install.packages("randomForest")
# Carregando o pacote
library(randomForest)
# Mensagem esperada:
# randomForest 4.7-1.1
# Type rfNews() to see new features/changes/bug fixes.

# Instalando o pacote
install.packages("fastICA")
# Carregando o pacote
library(fastICA)
# Mensagem esperada:
# (Sem mensagem)

# Instalando o pacote
install.packages("vegan")
# Carregando o pacote
library(vegan)
# Mensagem esperada:
# Carregando pacotes exigidos: permute
# Carregando pacotes exigidos: lattice
# This is vegan 2.6-2

# Instalando o pacote
install.packages("Rtsne")
# Carregando o pacote
library(Rtsne)
# Mensagem esperada:
# (Sem mensagem)

# Instalando o pacote
install.packages("umap")
# Carregando o pacote
library(umap)
# Mensagem esperada:
# (Sem mensagem)

# Instalando o pacote
install.packages('corrplot')
# Carregando o pacote
library(corrplot)
# Mensagem esperada:
# (Sem mensagem)

# Instalando o pacote
install.packages('factoextra')
# Carregando o pacote
library(factoextra)
# Mensagem esperada:
# (Sem mensagem)

# Instalando o pacote
install.packages("scatterplot3d")
# Carregando o pacote
library(scatterplot3d)
# Mensagem esperada:
# (Sem mensagem)
