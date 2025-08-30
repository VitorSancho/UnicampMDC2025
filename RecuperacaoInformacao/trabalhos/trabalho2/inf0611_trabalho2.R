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
class(imagens[[2]][2])

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
  lbp(grayscale(img)[,,1,1],1)
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

features_t <- t(sapply(imagens, lbp_desc))
rownames(features_t) <- names(imagens)

features_s <- t(sapply(imagens, Momentos))
rownames(features_s) <- names(imagens)

#----------------------------------------------------------------#
# Questao 2                               
#----------------------------------------------------------------#

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

# plotando as imagens retornadas
par(mfrow = c(3,3), mar = rep(2, 4))
plot_ranking <- function(query, features){
  distancia <- dist(features, method = "euclidean")
  distancia <- as.matrix(distancia)
  distancia_interesse <- distancia[,query]
  
  ranking <- order(distancia_interesse)
  return(ranking)
}
#-----------------------------#
# construindo rankings                          
# para cada uma das 5 consultas, construa um ranking com base na cor
ranking_c_biloba <- plot_ranking(consulta_biloba, features_c)
ranking_c_europaea <- plot_ranking(consulta_europaea, features_c)
ranking_c_ilex <- plot_ranking(consulta_ilex, features_c)
ranking_c_monogyna <- plot_ranking(consulta_monogyna, features_c)
ranking_c_regia <- plot_ranking(consulta_regia, features_c)

# para cada uma das 5 consultas, construa um ranking com base na textura
ranking_t_biloba <- plot_ranking(consulta_biloba, features_t)
ranking_t_europaea <- plot_ranking(consulta_europaea, features_t)
ranking_t_ilex <- plot_ranking(consulta_ilex, features_t)
ranking_t_monogyna <- plot_ranking(consulta_monogyna, features_t)
ranking_t_regia <- plot_ranking(consulta_regia, features_t)
  
# para cada uma das 5 consultas, construa um ranking com base na forma
ranking_s_biloba <- plot_ranking(consulta_biloba, features_s)
ranking_s_europaea <- plot_ranking(consulta_europaea, features_s)
ranking_s_ilex <- plot_ranking(consulta_ilex, features_s)
ranking_s_monogyna <- plot_ranking(consulta_monogyna, features_s)
ranking_s_regia <- plot_ranking(consulta_regia, features_s)

#-----------------------------#
# comparando  rankings                              
s
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
# (e) 
#                                         
#                                         
#                                         
#                                         
#                                         
#                                         
# (f)  
#                                         
#                                         
#                                         
#                                         
#                                         
# (g) dica: use a função generate_df_11_points para gerar o 
# data.frame interpolado para cada descritor. Depois, use a função
# plot_precision_x_recall_11_points_t2 para gerar o gráfico. 
# Exemplo de chamada da função generate_df_11_points:
# df_c_11_points <- generate_df_11_points(list(list(ground_truth_biloba, ranking_c_biloba), list(ground_truth_europaea, ranking_c_europaea)))
# Exemplo de chamada da função plot_precision_x_recall_11_points_t2:
# plot_precision_x_recall_11_points_t2(rbind(df_c_11_points, df_t_11_points, df_s_11_points), c("Cor", "Textura", "Forma"), "Titulo do Gráfico") 
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
ranking_concat_biloba   <- plot_ranking(consulta_biloba, features_concat)
ranking_concat_europaea <- plot_ranking(consulta_europaea, features_concat)
ranking_concat_ilex     <- plot_ranking(consulta_ilex, features_concat)
ranking_concat_monogyna <- plot_ranking(consulta_monogyna, features_concat)
ranking_concat_regia    <- plot_ranking(consulta_regia, features_concat)
  
# analisando rankings gerados com caracteristicas concatenadas
concat_biloba_analyse   <- analyse_rankings(ranking_concat_biloba, ground_truth_biloba)
concat_europaea_analyse <- analyse_rankings(ranking_concat_europaea, ground_truth_europaea)
concat_ilex_analyse     <- analyse_rankings(ranking_concat_ilex, ground_truth_ilex)
concat_monogyna_analyse <- analyse_rankings(ranking_concat_monogyna, ground_truth_monogyna)
concat_regia_analyse    <- analyse_rankings(ranking_concat_regia, ground_truth_regia)


#----------------------------------------------------------------#
# Questao 3 - RESPONDA:  
# (d) 
# 
# 
# (e) 
# 
# 
# 
# 
# (f)
# 
# 
# 
#----------------------------------------------------------------#




#----------------------------------------------------------------#
# Questao 4
#----------------------------------------------------------------#

# Definindo a consulta (mesmo índice da Questão 2)
consulta <- 4

# calculando as distancias, descritor:  histograma de cor 
dist_hist_q4 <- get_distance_vector(histogramas, consulta) 
r_hist_q4 <- order(dist_hist_q4)
  
# calculando as distancias, descritor:  textura 
dist_text_q4 <- get_distance_vector(texturas, consulta) 
r_text_q4 <- order(dist_text_q4)
  
# calculando as distancias, descritor:  forma 
dist_forma_q4 <- get_distance_vector(moments, consulta) 
r_forma_q4 <- order(dist_forma_q4)
  
# calculando e analisando rankings combmax
r_combmax_q4 <- names(imagens)[combmax(dist_hist_q4, dist_text_q4, dist_forma_q4)]
r_combmax_q4

# calculando e analisando rankings combsum
r_combsum_q4 <- names(imagens)[combsum(dist_hist_q4, dist_text_q4, dist_forma_q4)]
r_combsum_q4

# calculando e analisando rankings borda
r_borda_q4 <- names(imagens)[borda(r_hist_q4, r_text_q4, r_forma_q4)]
r_borda_q4

  
# analisar resultados
analyse_rankings(r_combmin_q4)
analyse_rankings(r_combmax_q4)
analyse_rankings(r_combsum_q4)
analyse_rankings(r_borda_q4)
  
#----------------------------------------------------------------#
# Questao 4 - RESPONDA:                   
# (i) 
# 
# 
# 
# (j)
# 
# 
# 
# (k)
# 
# 
#
# (l)
# 
#
#
#----------------------------------------------------------------#
