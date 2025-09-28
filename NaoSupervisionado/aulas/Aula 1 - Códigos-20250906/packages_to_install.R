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
# Pacotes para a Aula 1
##################################################################
# Instalando o pacote 
install.packages("caret")
# Carregando o pacote
library(caret)
# Mensagem esperada:
# Loading required package: lattice
# Loading required package: ggplot2

# Instalando o pacote 
install.packages("e1071")
# Carregando o pacote
library(e1071)
# Mensagem esperada:
# (Sem mensagem)

# Instalando o pacote 
install.packages("ggplot2")
# Carregando o pacote
library(ggplot2)
# Mensagem esperada:
# (Sem mensagem)

# Instalando o pacote 
install.packages("reshape2")
# Carregando o pacote
library(reshape2)
# Mensagem esperada:
# (Sem mensagem)


# Instalando o pacote 
install.packages("dplyr")
# Carregando o pacote
library(dplyr)
# Mensagem esperada:
# Attaching package: ‘dplyr’
#
# The following objects are masked from ‘package:stats’:
#  
#   filter, lag
#
# The following objects are masked from ‘package:base’:
#  
#   intersect, setdiff, setequal, union


# Instalando o pacote 
install.packages("MASS")
# Carregando o pacote
library(MASS)
# Mensagem esperada:
# Attaching package: ‘MASS’

# The following object is masked from ‘package:dplyr’:

#     select


# Instalando o pacote 
install.packages("mltools")
# Carregando o pacote
library(mltools)
# Possível mensagem esperada:
# Attaching package: ‘mltools’

# The following object is masked from ‘package:e1071’:

#     skewness


# Instalando o pacote 
install.packages("data.table")
# Carregando o pacote
library(data.table)
# Mensagem esperada:
# data.table 1.13.6 using 4 threads (see ?getDTthreads).  
# Latest news: r-datatable.com

# Attaching package: ‘data.table’

# The following objects are masked from ‘package:reshape2’:
  
#   dcast, melt

# The following objects are masked from ‘package:dplyr’:
  
#   between, first, last

# Instalando o pacote 
install.packages("nortest")
# Carregando o pacote
library(nortest)
# Mensagem esperada:
# (Sem mensagem)


# Instalando o pacote 
install.packages("outliers")
# Carregando o pacote
library(outliers)
# Mensagem esperada:
# (Sem mensagem)

# Instalando o pacote 
install.packages("EnvStats")
# Carregando o pacote
library(EnvStats)
# Mensagem esperada:
# 
# Attaching package: ‘EnvStats’

# The following object is masked from ‘package:mltools’:
  
#   skewness

# The following object is masked from ‘package:MASS’:
  
#   boxcox

# The following objects are masked from ‘package:stats’:
  
#   predict, predict.lm

# The following object is masked from ‘package:base’:
  
#   print.default


# Instalando o pacote 
install.packages("stringr")
# Carregando o pacote
library(stringr)
# Mensagem esperada:
# (Sem mensagem)




# Instalando o pacote 
install.packages("lubridate")
# Carregando o pacote
library(lubridate)
# Mensagem esperada:
# Attaching package: ‘lubridate’

# The following objects are masked from ‘package:data.table:
  
#   hour, isoweek, mday, minute, month, quarter, second, wday, week,
# yday, year

# The following objects are masked from ‘package:base’:
  
#   date, intersect, setdiff, union


# Instalando o pacote 
install.packages("tm")
# Carregando o pacote
library(tm)
# Mensagem esperada:
# Loading required package: NLP

# Attaching package: ‘NLP’

# The following object is masked from ‘package:ggplot2’:
  
#   annotate

# Instalando o pacote 
install.packages("wordcloud")
# Carregando o pacote
library(wordcloud)
# Mensagem esperada:
# Loading required package: RColorBrewer

# Instalando o pacote 
install.packages("textstem")
# Carregando o pacote
library(textstem)
# Mensagem esperada:
# Loading required package: koRpus.lang.en
# Loading required package: koRpus
# Loading required package: sylly
# For information on available language packages for 'koRpus', run

# available.koRpus.lang()

# and see ?install.koRpus.lang()