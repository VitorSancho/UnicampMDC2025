  #----------------------------------------------------------------#
# INF-0611 Recuperacao de Informacao       
#                       
# Trabalho Avaliativo 2
#----------------------------------------------------------------#
# Nome COMPLETO dos integrantes dp grupo:  
# - Vitor de Oliveira Fernandez Araujo                                       
# - Vitor Sancho Cardoso                                       
# -                                        
# 
#----------------------------------------------------------------#

#----------------------------------------------------------------#
# Configuracao dos arquivos auxiliares 
#----------------------------------------------------------------#
# configure o caminho antes de executar
# setwd("") 
options(warn=-1)
# install.packages("imager")
source("./ranking_metrics.R")
source("./trabalho2_base.R")

# caminho da pasta de imagens
path_plantas = "./plantas"

#----------------------------------------------------------------#
# Leitura das imagens 
#----------------------------------------------------------------#
imagens <- read_images(path_plantas)
class(imagens[[2]][0])


#----------------------------------------------------------------#
# Obtem classe de cada imagem 
#----------------------------------------------------------------#
nome_classes <- get_classes(path_plantas)
name_plantas <- list.files(path_plantas, full.names = FALSE)



#----------------------------------------------------------------#
# obtem ground_truth para cada classe 
#----------------------------------------------------------------#
ground_truth_biloba <- get_ground_truth(path_plantas, nome_classes, "biloba")
ground_truth_europaea <- get_ground_truth(path_plantas, nome_classes, "europaea")
ground_truth_ilex <- get_ground_truth(path_plantas, nome_classes, "ilex")
ground_truth_monogyna <- get_ground_truth(path_plantas, nome_classes, "monogyna")
ground_truth_regia <- get_ground_truth(path_plantas, nome_classes, "regia")


#----------------------------------------------------------------#
# Questao 1 
#----------------------------------------------------------------#

# obtem caracteristicas de cor  
hist_cor_desc <- function(path_img){
    img <- load.image(path_img)
    r <- hist(img[,,1]*255, plot=FALSE, breaks=0:255)$counts
    g <- hist(img[,,2]*255, plot=FALSE, breaks=0:255)$counts
    b <- hist(img[,,3]*255, plot=FALSE, breaks=0:255)$counts
    return(c(r, g, b))
}

# obtem caracteristicas de textura   
lbp_desc <- function(img){
  r1 <- lbp(grayscale(img)[,,1,1],1)
  lbp_uniforme <- hist(r1$lbp.u2, plot=FALSE, breaks=59)$counts
  return(c(lbp_uniforme))
}


# obtem caracteristicas de forma 
Momentos <- function(img){
  
  centroide <- function(M) {
    c(momento(M, 1, 0) / momento(M, 0, 0),
      momento(M, 0, 1) / momento(M, 0, 0))
  }
  
  momento <- function(M, p, q, central = FALSE) {
    r <- 0
    if (central) {
      c <- centroide(M)
      x <- c[1]
      y <- c[2]
    } else {
      x <- 0
      y <- 0
    }
    for (i in 1:nrow(M))
      for (j in 1:ncol(M))
        r <- r + (i - x)^p * (j - y)^q * M[i,j]  
    return(r)
  }
  
  img <- grayscale(img)[,,1,1]
  features <-NULL
  for(i in 0:2){
    for(j in 0:2){
      features <- cbind(features,momento(img, i,j, central=TRUE))
    }
  }
  return(features)
}


#----------------------------------------------------------------#
# obtem características de cor, textura e forma  
# para todas as imagens e armazena em matrizes 
# onde uma linha e uma imagem 
features_c <- t(sapply(names(imagens), hist_cor_desc))
rownames(features_c) <- names(imagens)

#features_t <- t(sapply(imagens, lbp_desc))
# pegando de csv para evitar tempo de processamento
features_t <- read.csv("features_t.csv")
rownames(features_t) <- names(imagens)

#features_s <- t(sapply(imagens, Momentos))
# pegando de csv para evitar tempo de processamento
features_s <- read.csv("features_s.csv")
rownames(features_s) <- names(imagens)

#----------------------------------------------------------------#
# Questao 2                               
#----------------------------------------------------------------#

## visualização das imagens retornadas por um ranking, para uma consulta
print_top_k <- function(query, ranking, k){
  
  par(mfrow = c(k+1,4), mar = rep(1, 4))
  plot(load.image(query), axes = FALSE, main = paste("consulta: ", query))
  for(img in ranking[1:k]){
    plot(load.image(img), axes = FALSE,
         main = img)
  }
}

# definindo as consultas
# obs.:  use o caminho completo para a imagem
consulta_biloba <- "./plantas/biloba_02.jpg"     
consulta_europaea <- "./plantas/europaea_01.jpg" 
consulta_ilex <- "./plantas/ilex_08.jpg"
consulta_monogyna <- "./plantas/monogyna_04.jpg"
consulta_regia <- "./plantas/regia_07.jpg"

# visualizando as consultas
par(mfrow = c(3,3), mar = rep(2, 4))
mostrarImagemColorida(consulta_biloba, "biloba_02.jpg")
mostrarImagemColorida(consulta_europaea,"europaea_01.jpg")
mostrarImagemColorida(consulta_ilex,"ilex_08.jpg")
mostrarImagemColorida(consulta_monogyna,"monogyna_04.jpg")
mostrarImagemColorida(consulta_regia,"regia_07.jpg")

#-----------------------------#
# construindo rankings                          
# para cada uma das 5 consultas, construa um ranking com base na cor
ranking_c_biloba <- get_ranking_by_distance(features_c, consulta_biloba)
ranking_c_europaea <- get_ranking_by_distance(features_c, consulta_europaea)
ranking_c_ilex <- get_ranking_by_distance(features_c, consulta_ilex)
ranking_c_monogyna <- get_ranking_by_distance(features_c, consulta_monogyna)
ranking_c_regia <- get_ranking_by_distance(features_c, consulta_regia)

# para cada uma das 5 consultas, construa um ranking com base na textura
ranking_t_biloba <- get_ranking_by_distance(features_t, consulta_biloba)
ranking_t_europaea <- get_ranking_by_distance(features_t, consulta_europaea)
ranking_t_ilex <- get_ranking_by_distance(features_t, consulta_ilex)
ranking_t_monogyna <- get_ranking_by_distance(features_t, consulta_monogyna)
ranking_t_regia <- get_ranking_by_distance(features_t, consulta_regia)
  
# para cada uma das 5 consultas, construa um ranking com base na forma
ranking_s_biloba <- get_ranking_by_distance(features_s, consulta_biloba)
ranking_s_europaea <- get_ranking_by_distance(features_s, consulta_europaea)
ranking_s_ilex <- get_ranking_by_distance(features_s, consulta_ilex)
ranking_s_monogyna <- get_ranking_by_distance(features_s, consulta_monogyna)
ranking_s_regia <- get_ranking_by_distance(features_s, consulta_regia)

#-----------------------------#
# comparando  rankings                              

## utilize as funções do arquivo ranking_metrics.R para calcular 
# a precisão, revocação, taxa F1 e precisão média nos 
# top 5, 10, 15 e 20
  
analyse_rankings <- function(ranking, ground_truth) {
  column_names <- c("top","Precisao", "Recall", "F1","AP")
  
  df <- data.frame(matrix(nrow = 0, ncol = length(column_names)))
  colnames(df) <- column_names
  
  for (x in c(5,10,15,20)) {
   df <- rbind(df, data.frame(top = paste("top_",as.character(x), sep=""),
                              Precisao = precision(ground_truth, ranking, x),
                              Recall = recall(ground_truth, ranking, x),
                              F1 = f1_score(ground_truth, ranking, x),
                              AP = ap(ground_truth, ranking, x)))
  }
  return(df)
}


# analisando rankings gerados com caracteristicas de cor
c_biloba_analyse <- analyse_rankings(ranking_c_biloba, ground_truth_biloba)
c_europaea_analyse <- analyse_rankings(ranking_c_europaea, ground_truth_europaea)
c_ilex_analyse <- analyse_rankings(ranking_c_ilex, ground_truth_ilex)
c_monogyna_analyse <- analyse_rankings(ranking_c_monogyna, ground_truth_monogyna)
c_regia_analyse <- analyse_rankings(ranking_c_regia, ground_truth_regia)

# analisando rankings gerados com caracteristicas de textura
t_biloba_analyse <- analyse_rankings(ranking_t_biloba, ground_truth_biloba)
t_europaea_analyse <- analyse_rankings(ranking_t_europaea, ground_truth_europaea)
t_ilex_analyse <- analyse_rankings(ranking_t_ilex, ground_truth_ilex)
t_monogyna_analyse <- analyse_rankings(ranking_t_monogyna, ground_truth_monogyna)
t_regia_analyse <- analyse_rankings(ranking_t_regia, ground_truth_regia)

# analisando rankings gerados com caracteristicas de forma
s_biloba_analyse <- analyse_rankings(ranking_s_biloba, ground_truth_biloba)
s_europaea_analyse <- analyse_rankings(ranking_s_europaea, ground_truth_europaea)
s_ilex_analyse <- analyse_rankings(ranking_s_ilex, ground_truth_ilex)
s_monogyna_analyse <- analyse_rankings(ranking_s_monogyna, ground_truth_monogyna)
s_regia_analyse <- analyse_rankings(ranking_s_regia, ground_truth_regia)

#----------------------------------------------------------------#
# Questao 2 - RESPONDA:                   
#############################################

## IMPORTANTE: para printar corretamente, precisa aumentar o tamanho da janela de plot 
print_top_k(consulta_europaea, ranking_c_europaea, 10)
print_top_k(consulta_europaea, ranking_t_europaea, 10)
print_top_k(consulta_europaea, ranking_s_europaea, 10)

c_europaea_analyse
t_europaea_analyse
s_europaea_analyse

#testando para encontrar K onde descritor Textura entrega Recall == 1
recall(ground_truth_europaea, ranking_t_europaea, 12)
recall(ground_truth_europaea, ranking_t_europaea, 13)

#############################################
# (e) 
#  Ao analisar a consulta com a imagem europaea_01.jpg, podemos perceber visualmente
#  que o descritor que performa melhor é o de forma. Intuitivamente, podemos perceber
#  que este descritor tem maior facilidade pois a folha tem um formato alongado, bem 
#  característico, dentre as 5 espécies analisadas. Ao analisar as métricas de avaliação
#  para os 3 descritores, o resultado matemático corrobora com a intuição, com o descritor
#  de forma sendo o único a conseguir precisão média 1 para todos os top-K, além de entregar
#  Recall=1 desde o top-10, feito que os descritores de cor e textura só conseguiram nos tops
#  20 e 13, respectivamente.
#                                         
# ###########################################

#calculando MAP para os descritores
color_map_list <- list(
  list(ground_truth_biloba, ranking_c_biloba),
  list(ground_truth_europaea, ranking_c_europaea),
  list(ground_truth_ilex, ranking_c_ilex),
  list(ground_truth_monogyna, ranking_c_monogyna),
  list(ground_truth_regia, ranking_c_regia)
)
texture_map_list <- list(
  list(ground_truth_biloba, ranking_t_biloba),
  list(ground_truth_europaea, ranking_t_europaea),
  list(ground_truth_ilex, ranking_t_ilex),
  list(ground_truth_monogyna, ranking_t_monogyna),
  list(ground_truth_regia, ranking_t_regia)
)

shape_map_list <- list(
  list(ground_truth_biloba, ranking_s_biloba),
  list(ground_truth_europaea, ranking_s_europaea),
  list(ground_truth_ilex, ranking_s_ilex),
  list(ground_truth_monogyna, ranking_s_monogyna),
  list(ground_truth_regia, ranking_s_regia)
)

MAP_color <- map(color_map_list, 10) ## 0.669
MAP_texture <- map(texture_map_list, 10) ## 0.492
MAP_shape <- map(shape_map_list, 10) ## 0.391

########################################
# (f) 
#  Com base na métrica de média das precisões médias em top 10, o descritor de cor obteve
#  melhor performance, atingindo uma MAP de 0.66. Para as 5 consultas, este descritor
#  foi o que teve a maior capacidade, em média, de retornar resultados relevantes no topo do ranking.
#                                         
######################################
#


df_c_11_points <- generate_df_11_points(
  list(
    list(ground_truth_biloba, ranking_c_biloba), 
    list(ground_truth_europaea, ranking_c_europaea),
    list(ground_truth_ilex, ranking_c_ilex),
    list(ground_truth_monogyna, ranking_c_monogyna),
    list(ground_truth_regia, ranking_c_regia)
  )
)
df_t_11_points <- generate_df_11_points(
  list(
    list(ground_truth_biloba, ranking_t_biloba),
    list(ground_truth_europaea, ranking_t_europaea),
    list(ground_truth_ilex, ranking_t_ilex),
    list(ground_truth_monogyna, ranking_t_monogyna),
    list(ground_truth_regia, ranking_t_regia)
  )
)
df_s_11_points <- generate_df_11_points(
  list(
    list(ground_truth_biloba, ranking_s_biloba),
    list(ground_truth_europaea, ranking_s_europaea),
    list(ground_truth_ilex, ranking_s_ilex),
    list(ground_truth_monogyna, ranking_s_monogyna),
    list(ground_truth_regia, ranking_s_regia)
  )
)

plot_precision_x_recall_11_points_t2(rbind(df_c_11_points, df_t_11_points, df_s_11_points), c("Cor", "Textura", "Forma"), "Curva PR - Precisão Média Interpolada em 11 pontos")
#
#
# (g) 
#  Analisando as curvas geradas, o descritor de cor se mostra muito melhor que os outros,
#  para niveis baixos de revocação. Conforme chegamos em niveis de revocação mais altos,
#  o descritor de textura atinge uma precisão média bem próxima, mas cai ao atingir revocação = 1.
#  Ao final, o descritor de cor consegue uma precisão média de 0.65, enquanto os descritores
#  textura e forma atingem 0.51 e 0.41, respectivamente, mostrando uma maior capacidade da cor
#  em caracterizar os elementos do conjunto em consultas para múltiplas espécies. Este resultado
#  corrobora com a quantidade de features em cada descritor, dado que o descritor de forma tem
#  apenas 9 features, o de textura (LBP uniforme) tem 59 e o de cor tem 765 (255*3). Isto faz com que
#  a característica de cor consiga codificar mais nuances do conjunto de dados, trazendo também um
#  resultado agregado mais preciso.
#                                         
#                                         
#                                         
#----------------------------------------------------------------#



#----------------------------------------------------------------#
# Questao 3
#----------------------------------------------------------------#
# concatenando caracteristicas                      

## obter vetores finais de caracteristicas pela concatenação de 
# cada tipo de caracteristica (cor, textura e forma):
features_concat <- cbind(features_c, features_t, features_s)
  
# gerar novos rankings
ranking_concat_biloba   <- get_ranking_by_distance(features_concat, consulta_biloba)
ranking_concat_europaea <- get_ranking_by_distance(features_concat, consulta_europaea)
ranking_concat_ilex     <- get_ranking_by_distance(features_concat, consulta_ilex)
ranking_concat_monogyna <- get_ranking_by_distance(features_concat, consulta_monogyna)
ranking_concat_regia    <- get_ranking_by_distance(features_concat, consulta_regia)
  
# analisando rankings gerados com caracteristicas concatenadas
concat_biloba_analyse   <- analyse_rankings(ranking_concat_biloba, ground_truth_biloba)
concat_europaea_analyse <- analyse_rankings(ranking_concat_europaea, ground_truth_europaea)
concat_ilex_analyse     <- analyse_rankings(ranking_concat_ilex, ground_truth_ilex)
concat_monogyna_analyse <- analyse_rankings(ranking_concat_monogyna, ground_truth_monogyna)
concat_regia_analyse    <- analyse_rankings(ranking_concat_regia, ground_truth_regia)


#----------------------------------------------------------------#
# Questao 3 - RESPONDA:  
# (d) 
library(tidyverse)
s_regia_analyse_para_comparacao <- s_regia_analyse 
names(s_regia_analyse_para_comparacao) <- c("top ", "Precisao_s", "Recall_s", "F1_s", "AP_s")
t_regia_analyse_para_comparacao <- t_regia_analyse
names(t_regia_analyse_para_comparacao) <- c("top ","Precisao_t", "Recall_t", "F1_t", "AP_t")
c_regia_analyse_para_comparacao <- c_regia_analyse
names(c_regia_analyse_para_comparacao) <- c("top ","Precisao_c", "Recall_c", "F1_c", "AP_c")
concat_regia_analyse

# Utilizamos as consultas de regia para realizar a comparação
# de desempenho entre os descritores isolados e o combinado.
# Comparando os resultados observa-se que a combinação entre os descritores
# teve uma piora no resultado se comparado com as features individuais.
# A tendência já observada se manteve, que o descritor de textura foi
# o mais eficiente, seguido do de cor e por último o combinado empatado com
# o de forma. O descritor combinado acabou se igualando ao descritor 
# menos informativo, em termos isolados.

# (e) 

# Ao analisar o resultado do descritor combinado é possível notar que
# seu resultado está igual ao obtido para para o descritor de forma. O que
# indica que este descritor é muito mais influente na qualidade do resultado
# do que os demais, tendo um maior poder informativo ao se utilizar features
# combinadas. Esse resultado se justifica ao olhar as imagens de cada
# uma das espécies. Estas possuem formas muito distintas entre si. Esse 
# resultado em um primeiro momento parece estranho, pois isoladamente, foi
# o descritor de pior resulado. Isso indica que de forma isolada este não é
# capaz de causar diferenciação entre as amostras, mas ao se combinar
# multiplas features, este foi o que teve mais poder de diferenciação
# entre as amostrs.

# (f)
# 
# 
# 
#----------------------------------------------------------------#




#----------------------------------------------------------------#
# Questao 4
#----------------------------------------------------------------#

# Definindo a consulta (mesmo índice da Questão 2)
consulta_regia

# calculando as distancias, descritor:  histograma de cor 
dist_hist_q4 <- get_distance_vector(features_c, consulta_regia) 
r_hist_q4 <- order(dist_hist_q4)
  
# calculando as distancias, descritor:  textura 
dist_text_q4 <- get_distance_vector(features_t, consulta_regia) 
r_text_q4 <- order(dist_text_q4)
  
# calculando as distancias, descritor:  forma 
dist_forma_q4 <- get_distance_vector(features_s, consulta_regia) 
r_forma_q4 <- order(dist_forma_q4)
  
# calculando e analisando rankings combmax
r_combmax_q4 <- names(imagens)[combmax(dist_hist_q4, dist_text_q4, dist_forma_q4)]
r_combmax_q4

# calculando e analisando rankings combsum
r_combmin_q4 <- names(imagens)[combmin(dist_hist_q4, dist_text_q4, dist_forma_q4)]
r_combmin_q4

# calculando e analisando rankings combsum
r_combsum_q4 <- names(imagens)[combsum(dist_hist_q4, dist_text_q4, dist_forma_q4)]
r_combsum_q4

# calculando e analisando rankings borda
r_borda_q4 <- names(imagens)[bordacount(r_hist_q4, r_text_q4, r_forma_q4)]
r_borda_q4

  
# analisar resultados
analyse_rankings(r_combmax_q4, ground_truth_regia)
analyse_rankings(r_combmin_q4, ground_truth_regia)
analyse_rankings(r_combsum_q4, ground_truth_regia)
analyse_rankings(r_borda_q4, ground_truth_regia)
  
#----------------------------------------------------------------#
# Questao 4 - RESPONDA:                   
# (i) 
# Para a consulta selecionada (regia), o método combsum retornou o melhor ranking.
# Observando os resultados é possível notar um bom equilíbro entre precisão e o
# recall. Para o top 5 e top 10 o desempenho é melhor em comparação com os outros,
# tendo uma precisão de 0.8 e 0.7 e um recall de 0.4 e 0.7, respectivamente. O método
# combMax teve uma performance muito parecida ao comparar os resultados para os top 15
# e 20 os resultados, mas os de top 5 e 10 demontram como esse método não foi capaz 
# de retornar melhores resultados, sendo o combSum o que retornou itens corretor de forma "mais rápida".
# 
# (j)
# 
consultas <- list(
list(ground_truth_regia, list(r_combmax_q4, r_combmin_q4, r_combsum_q4, r_borda_q4))
)
df_combmax <- generate_df_11_points(list(list(ground_truth_regia, r_combmax_q4)))
df_combmin <- generate_df_11_points(list(list(ground_truth_regia, r_combmin_q4)))
df_combsum <- generate_df_11_points(list(list(ground_truth_regia, r_combsum_q4)))
df_borda   <- generate_df_11_points(list(list(ground_truth_regia, r_borda_q4)))

df_metodos <- rbind(
  df_combmax, df_combmin, df_combsum, df_borda
)

metodos <- c("CombMAX", "CombMIN", "CombSUM", "BORDA")

plot_precision_x_recall_11_points_t2(df_metodos, metodos, "Curva Precisão x Revocação - Métodos de Agregação")

# Ao observar o gráfico obtido, é possível confirmar o que foi constato na questão anterior,
# sendo o método CombSum o de melhor desempenho. Porém agora é possível notar a similaridade
# entre o CombSum eo CombMax. Ambos só diferem em precisão para uma revocação de 0.6. Para
# os demais pontos os métodos são equivalentes.


# (k)
# 
# 
#
# (l)
# 
#
#
#----------------------------------------------------------------#
