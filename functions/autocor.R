# ------------------------------------------------------------ #
# Rotina para geracao de mapas tematicos, moran scatterplots,
# box maps, lisa maps e Moran maps e calculo dos indices global
# e local de Moran.
#
# Autor: Leonardo Gravina de Faria
#
# Desenvolvido como projeto de iniciacao cientifica na Univer-
# sidade Federal de Juiz de Fora entre agosto de 2018 e julho
# de 2019 sob orientacao do professor Tiago Maia Magalhaes.
#
# Os codigos sao baseados naqueles utilizados pelo professor 
# na disciplina de Introducao a Estatistica Espacial em 2018-2.
# 
# Ultima versao: 16/10/2019
# ------------------------------------------------------------ #

library(dplyr)
library(rgdal)
library(ggplot2)
library(spdep)
library(Cairo)

# shp: variavel contendo o objeto do tipo SpatialPolygonsDataFrame
# id.varb: atributo de identificacao individual de cada poligono
# data: uma lista com os dados que deve conter os seguintes itens:
    # varb: a variavel analisada
    # varname: o nome da variavel com que os graficos serao salvos. 
    #  Ex: "TEM_varname.png" para o mapa tematico.
    # vartitle: o nome da variavel como titulo do grafico
    # varcolor: a palheta de cores continua com que o mapa tematico sera preenchido
    # interval: vetor com os limites dos intervalos para o mapa tematico
# matrix.mode: como a matriz de vizinhanca sera preenchida na funcao nb2mat. 
#  "W" para matriz de pesos, "B" para matriz binaria.
# sign.lvl: nivel de significancia com que sera gerado o Moran map
# width e heigth: dimensoes da imagem gerada
# ggadd: objeto ggplot adcicional para ser adicionado ao plot
#
# retorno: lista com o teste de I global de Moran e I locais

moran.autocor = function(shp, id.varb = 'id', data, matrix.mode = 'W', mNB = NULL, sign.lvl = .05, width = 5, height = 4.5, ggadd = NULL){
  # ggplot2 gera mapas mais facilmente se o SpatialPolygonsDataFrame for
  # convertido em um dataframe. A funcao fortify faz isso.
  mTract = fortify(shp, region = id.varb)
  mTract$id = as.character(mTract$id)
  if(is.null(mNB))
    mNB = poly2nb(shp) #gera uma lista de vizinhancas
  mW = nb2mat(mNB, style = matrix.mode) #gera a matriz de vizinhancas
  
  cat(data$varname,':\nMatriz de vizinhanca montada\n')
  
  varname = data$varname
  vartitle = data$vartitle
  varcolor = data$varcolor
  varb = data$varb
  
  id = unique(mTract$id)
  
  # ----- Mapa Tematico ----- #
  
  vCuts = cut(varb, data$interval, includes.lowest = TRUE, ordered_result = TRUE) 
  vLevels = levels(vCuts)
  
  # Une as areas do dataframe com os fatores correpondentes
  plotData = left_join(mTract, data.frame(cbind(id = as.character(id), IC = as.character(vCuts))))
  plotData$IC = factor(plotData$IC, levels = vLevels, ordered = T)
  
  cat('Mapa tematico\n')
  
  pTem = ggplot() +
    geom_polygon(data = plotData, aes(x = long, y = lat, group = group, fill = IC), color = gray(0.1), size = 0.05) +
    coord_map() +
    scale_fill_brewer(palette = varcolor, direction = 1) + ggadd +
    labs(title = vartitle, fill = "", x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  cat('Mapa tematico pronto\n')
  
  ggsave(pTem, file = paste0("./Outputs/TEM_", varname, ".png"), width = width, height = height, type = "cairo-png")
  
  cat('Mapa tematico salvo\n')
  # ----- Moran Scatter plot ----- #
  
  vZ = varb - mean(varb) # residuos
  vWZ = mW %*% vZ # soma dos residuos dos vizinhos
  
  # associa cada area ao respectivo quadrante
  vColor = NA
  vColor[which(vZ >= 0 & vWZ >= 0)] = "Q1"
  vColor[which(vZ <  0 & vWZ <  0)] = "Q2"
  vColor[which(vZ <  0 & vWZ >= 0)] = "Q3"
  vColor[which(vZ >= 0 & vWZ <  0)] = "Q4"
  
  mMoran = data.frame(vZ, vWZ, vColor)
  
  cat('Moran scatter plot\n')
  
  pMSP = ggplot(mMoran, aes(x = vZ, y = vWZ, colour = vColor)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    labs(title = vartitle, x = "Z", y = "WZ") +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank())
  
  cat('Moran scatter plot pronto\n')
  
  ggsave(pMSP, file = paste0("./Outputs/MoranSP_", varname, ".png"), width = width, height = height, type = "cairo-png")
  
  cat('Moran scatter plot salvo\n')  
  # ---- Box map
  
  plotDataBM = left_join(mTract, data.frame(cbind(id = as.character(id), Quad = vColor)))
  
  cat('plotDataBM\n')
  
  pBM = ggplot() +
    geom_polygon(data = plotDataBM, aes(x = long, y = lat, group = group, fill = Quad), color = gray(0.1), size = 0.05) +
    coord_map() +
    scale_fill_manual(breaks = c("Q1", "Q2", "Q3", "Q4"), 
                      values = c("firebrick", "dodgerblue2", "gold", "springgreen")) + ggadd +
    labs(title = vartitle, fill = "", x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  cat('Box map pronto\n')
  
  ggsave(pBM, file = paste0("./Outputs/Boxmap_", varname, ".png"), width = width, height = height, type = "cairo-png")
  
  cat('Box map salvo\n')
  
  # ---- Global Moran's I
  # Calculo sem uso da funcao
  # dS0 = sum(mW)
  # dI = (length(vZ) * t(vZ) %*% mW %*% vZ)/(dS0 * t(vZ) %*% vZ)
  # 
  # cat('Moran\'s Global I: ', dI, '\n')

  mWlist = nb2listw(mNB, style = matrix.mode)
  globalMoran = moran.test(vZ, mWlist)
  
  # ---- Local Moran's I
  # Calculo sem o uso da funcao
  # mZd = diag(x = 1, nrow = length(vZ))
  # diag(mZd) = vZ
  # vI = length(vZ) * (mZd %*% mW %*% vZ)/as.numeric(t(vZ) %*% vZ)
  # 
  # cat('Moran\'s Local I\'s:', vI, '\n')
  
  mLMoran = localmoran(vZ, mWlist)
  
  # ---- Lisa map
  
  # associa cada area a respectiva significancia do indice local de autocorrelacao espacial
  vColorLisa = NA
  vColorLisa[which(mLMoran[, 5] <= 0.001)] = "99.9%"
  vColorLisa[which(mLMoran[, 5] > 0.001 & mLMoran[, 5] <= 0.010)] = "99%"
  vColorLisa[which(mLMoran[, 5] > 0.010 & mLMoran[, 5] <= 0.050)] = "95%"
  vColorLisa[which(mLMoran[, 5] > 0.050)] = "NS"
  
  plotDataLM = left_join(mTract, data.frame(cbind(id = as.character(id), Quad = vColorLisa)))
  
  cat('LISA map\n')
  
  pLM = ggplot() +
    geom_polygon(data = plotDataLM, aes(x = long, y = lat, group = group, fill = Quad), color = "black", size = 0.05) +
    coord_map() +
    scale_fill_manual(breaks = c("NS", "95%", "99%", "99.9%"), 
                      values = c("NS" = "#FFFFFF", "95%" = "#FFFF00", "99%" = "#FF9900", "99.9%" = "#FF0000")) + ggadd +
    labs(title = vartitle, fill = "", x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  cat('LISA map pronto\n')
  
  ggsave(pLM, file = paste0("./Outputs/Lisamap_", varname, ".png"), width = width, height = height, type = "cairo-png")
  
  cat('LISA map salvo\n')
  
  # ---- Moran map
  
  vColorMoran = rep("NS", length(vColor))
  vColorMoran[which(vColor == "Q1" & mLMoran[, 5] <= sign.lvl)] = "Q1"
  vColorMoran[which(vColor == "Q2" & mLMoran[, 5] <= sign.lvl)] = "Q2"
  vColorMoran[which(vColor == "Q3" & mLMoran[, 5] <= sign.lvl)] = "Q3"
  vColorMoran[which(vColor == "Q4" & mLMoran[, 5] <= sign.lvl)] = "Q4"
  
  plotDataMM = left_join(mTract, data.frame(cbind(id = as.character(id), Quad = vColorMoran)))
  
  cat('Moran map\n')
  
  pMM = ggplot() +
    geom_polygon(data = plotDataMM, aes(x = long, y = lat, group = group, fill = Quad), color = "black", size = 0.05) +
    coord_map() +
    scale_fill_manual(breaks = c("NS", "Q1", "Q2", "Q3", "Q4"), 
                      values = c("white", "firebrick", "dodgerblue2", "gold", "springgreen")) + ggadd +
    labs(title = vartitle, fill = "", x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  cat('Moran map pronto\n')
  
  ggsave(pMM , file = paste0("./Outputs/Moranmap_", varname, ".png"), width = width, height = height, type = "cairo-png")
  
  cat('Moran map salvo\n')
  
  return(list(globalMoran = globalMoran, localMoran = mLMoran))
}

## ---- Funcoes auxiliares ---- ##

# Retorna uma lista no formato usado pela funcao principal para o dado analisado
var.info = function(varb = NULL, varname = "", vartitle = "", varcolor = "", interval = c(-Inf,+Inf)){
  if(class(interval)=='function')
    interval = interval(varb)
  return(list(varb = varb, varname = varname, vartitle = vartitle, varcolor = varcolor, interval = interval))
}

# Retorna um vetor com os limites dos intervalos
# iguais aos valores do
# minimo, limite inferior, primeiro quartil, mediana,
# terceiro quartil, limite superior e maximo.
quantileInt = function(vrb){
  vec = unname(quantile(vrb))
  limInf = vec[2] - 1.5*(vec[4]-vec[2])
  limSup = vec[4] + 1.5*(vec[4]-vec[2])
  retorno = c(vec[1] - .01)
  if(limInf > vec[1]){
    retorno[2:5] = c(limInf, vec[2:4])
    if(limSup < vec[5])
      retorno = c(retorno, limSup, vec[5])
    else
      retorno = c(retorno, vec[5])
  }
  else{
    if(limSup < vec[5])
      retorno[2:6] = c(vec[2:4], limSup, vec[5])
    else
      return(c(retorno, vec[2:5]))
  }
  return(retorno)
}
