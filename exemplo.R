# ------------------------------------------------------------ #
# Exemplo do uso da funcao para indices de moran e mapas 
# relativos. Os dados usados sao os socioeconomicos relativos
# ao censo 2010 para os municipios do estado de Minas Gerais.
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
# Caso se deseje realizar analises por micro ou mesorregiao, 
# o arquivo dtb_2010(MG).csv fornece os nomes e codigos das 
# mesmas.
# 
# Ultima versao: 22/10/2019
# ------------------------------------------------------------ #

source("functions/autocor.r")

MGSHP = readOGR(dsn = "areas" , layer = "31MUE250GC_SIR")
MGSHP = MGSHP[order(MGSHP$CD_GEOCODM),] # ordenar a malha pelo codigo do municipio

mNB = poly2nb(MGSHP) # lista de vizinhancas

## ----- Importando as malhas digitais ----- ##

# importando o shp para as mesorregioes
MG_MESO = readOGR(dsn = "areas", layer = "31MEE250GC_SIR")
MG_MESO = MG_MESO[order(MG_MESO$CD_GEOCODU),] # ordenar a malha pelo codigo das mesorregioes
mTract_meso = fortify(MG_MESO, region = "ID")

# Objeto ggplot para a malha de mesorregioes
plot_meso = geom_polygon(data = mTract_meso, 
                         aes(x = long, y = lat, group = group, fill = NA), 
                         color = "black", size = 0.1)

## ----- Importando os dados ----- ##

df = read.csv2("data/AtlasBrasil_IDHM_MG.csv", header = T)
df = df[order(df$id),] # ordenar os dados pelos codigos dos municipios

# Calculando a taxa de frequencia de matricula

df$frequencia = (df$freq.5a6 + 
  df$freq.11a13.fundamental.final + 
  df$freq.18a20.medio.completo + 
  df$freq.15a17.fundamental.completo) / 4

# --------------------------------------#
# Intervalos de IDH separados em 
# IDHM  < 0.500: Muito baixo
# 0.500 - 0.599: Baixo
# 0.600 - 0.699: Medio
# 0.700 - 0.799: Elevado
# 0.800 -     1: Muito elevado
#
idhmInt <- c(0, 0.499, 0.599, 0.699, 0.799, 1)
# --------------------------------------#

variaveis = list(IDHM = var.info(df$idhm,'idhm', 'IDHM', 'Blues', idhmInt),
                 RendaPC = var.info(df$renda, 'renda', 'Renda per capita', 'Greens', quantileInt),
                 ExpVida = var.info(df$exp.vida, 'expvida', 'Expectativa de vida em anos', 'Reds', quantileInt),
                 IDHM.edu = var.info(df$idhm.edu, 'edu', 'IDHM Educacional', 'Purples', idhmInt),
                 Escolaridade = var.info(df$escolaridade, 'escolaridade', 'Escolaridade', 'Purples', quantileInt),
                 Frequencia = var.info(df$frequencia, 'frequencia', 'Taxa de frequencia de matricula', 'Purples', quantileInt))

resultados = lapply(variaveis, function(X){
 moran.autocor(MGSHP, 'CD_GEOCODM', X, ggadd = plot_meso, width = 6, mNB = mNB)
 })

moran.autocor(MGSHP, 'CD_GEOCODM', variaveis$IDHM, ggadd = plot_meso, width = 6, mNB = mNB)

