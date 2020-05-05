
# limpando todas as variaveis antes de comecar
rm(list = ls())

#--------------------------------------------------------------------
# Estimacao de parametros - Modelo estocastico
#
#       beta           gamma
# S + I ----> 2I  |  I -----> R
#
# I = I_s + I_a,  onde I_a = m * I_s  
# - fixando m = 3
# - fixando gamma = 1
#
#--------------------------------------------------------------------

# PARAMETROS PARA A MONTAGEM DA BASE DE DADOS
qAno <- 40         # semana de quebra de ano
anoRef <- 2013     # ano de referencia

# PARAMETROS DO MODELO
m <- 3             # parametros que relaciona os assintomaticos
nSemAjuste <- 5    # numero de semanas de ajuste
pSuscPop <- 1.0    # proporcao de suscetiveis na populacao inicial
nSemFut <- 4       # numero de semanas no futuro
nMedMovel <- 3     # numero de semanas para suavizacao por medias moveis
tamMH <- 100       # tamanho da cadeia Metropolis-Hastings

# parametros do MCMC
# priori
# parametros da distribuicao gamma a priori
a1 <- 1.01
a2 <- 1.01
b1 <- 1
b2 <- 1

############################
# caminhos de bases/saidas #
############################
# caminho da base de notificacao
caminhoNot <- "/home/nakashima/Documentos/IME/Epidemio/R/extracaoParametros/notificacoes.csv"
# caminho da base de Populacao
caminhoPop <- "/home/nakashima/Documentos/IME/Epidemio/R/extracaoParametros/populacao.csv"
# caminho da base de nomes das cidades
caminhoNom <- "/home/nakashima/Documentos/IME/Epidemio/R/extracaoParametros/nomesCidades.txt"
# caminho para salvar grafico com semana de ajuste
caminhoAjuste <- "/home/nakashima/perAjusteB44.pdf"

#############################
# MONTAGEM DA BASE DE DADOS #
#############################
# inicio e final do ano no formato AAAASS
iniAno <- (anoRef - 1) * 100 + qAno
fimAno <- anoRef * 100 + qAno - 1

# notificacoes
Not <- read.csv(caminhoNot, header = FALSE, sep = " ")
colnames(Not) <- Not[1,]         # nome das colunas eh a primeira linha
Not <- Not[-1,]                  # eliminando a linha com o nome das colunas

indAno <- match(c(iniAno, fimAno), names(Not))     # encontrando inicio e fim do ano estudado
Not <- Not[, indAno[1]:indAno[2]]                  # olhando apenas para as notificacoes do ano em questao
Not <- t(Not)                                      # transpondo (facilita a vida)
colnames(Not) <- c(1:length(Not[1,]))              # id da cidade

Not_mAux <- Not * (1 + m)                          # infectados = sintomaticos + assintomaticos
Not_m <- Not_mAux                                  # infectados = (1 + m) * sintomaticos

nSemAno <- length(Not_m[,1])                       # numero de semanas do ano escolhido

# suavizacao por medias moveis
for(i in 3:nSemAno)
{
  Not_m[i, ] <- round((Not_mAux[i-2,] + Not_mAux[i-1,] + Not_mAux[i, ])/3, digits = 0)
}

# populacao
Pop <-read.csv(caminhoPop, header = FALSE, sep = " ")
Pop <- Pop[, 2]

# nomes
Nomes <- read.csv(caminhoNom, header = FALSE)
Nomes <- t(Nomes)                                  # transposto para facilitar

###########
#         #
# FUNCOES #
#         #
###########
#--------------------------------------------------------
# organiza os dados semanalmente no vetor de estados X(t)
# carregamento de X(t)
# X(t) = (S(t), I(t), R(t))
# X(t) vetor de estados para cada semana t e para cada cidade
# recebe: Not_m[semana, cidade] matriz de notificacoes * (1 + m) por semana e cidade
#         Pop[cid] populacao por cidade, nSemAjuste numero de semanas do ajuste, pSuscPop
#         in [0,1] proporcao de suscetiveis na populacao inicial, ind_t0[cid] vetor de
#         indices para semana inicial por cidade, cid[cid] vetor com as cidades estudadas
# devolve: X$?(S, I ou R) estado observado em tempo discreto
geraDados <- function(Not_m, Pop, nSemAjuste, pSuscPop, ind_t0, cid)
{
  S <- matrix(nrow = nSemAjuste, ncol = length(cid)) # suscetiveis[semana, cid]
  I <- matrix(nrow = nSemAjuste, ncol = length(cid)) # infectados[semana, cid]
  R <- matrix(nrow = nSemAjuste, ncol = length(cid)) # recuperados[semana, cid]
  
  for(i in 1:length(cid))
  {
    # primeira semana do periodo de ajuste
    S[1, i] <- round(Pop[cid[i]] * pSuscPop) - Not_m[ind_t0[i], i]
    I[1, i] <- Not_m[ind_t0[i], i]
    R[1, i] <- Pop[cid[i]] - round(Pop[cid[i]] * pSuscPop)
    
    # demais semanas do periodo de ajuste
    for(j in 2:nSemAjuste)
    {
      S[j, i] <- S[j-1, i] - Not_m[ind_t0[i]+j-1, i]
      I[j, i] <- Not_m[ind_t0[i]+j-1, i]
      R[j, i] <- R[j-1, i] + Not_m[ind_t0[i]+j-2, i]
    }
  }
  
  X <- list(S = S, I = I, R = R)
  return(X)
}

#-----------------------------------------------------------------------------------
# transformacao de tempo que relaciona os processos homogeneo e nao homogeneo linear
# recebe: tempo em processo de Poisson homogeneo
# devolve: tempo transformado - Processo de Poisson nao homogeneo linear
# recebe: h0, h1 valor das propensidades nos extremos de um intervalo, t[n] vetor de tempos
#         de ocorrencia de eventos eno intervalo, i indice do intervalo analizado
# devolve: tempo[n] tempo transformado
transformaTempo <- function(h0, h1, t, i)
{
  if(h0 == h1)
    return(t)
  else
    return(i + (sqrt(h0^2 + (h1^2 - h0^2) * (t - i)) - h0)/(h1 - h0))
}

#-----------------------------------------------------------------------------------
# gera os instantes de reacoes distribuindo os n eventos dentro de seu respectivo
# intervalo segundo a distribuição uniforme
# apos esse passo, utiliza a transformacao definida na funcao transformaTempo para
# transforma-lo em tempos provenientes de um processo de Poisson não uniforme
# recebe: S[sem, cid] e I[sem, cid] matrizes de suscetiveis por semana e cidade, cid[cid]
#         vetor com cidades estudadas, nSemAjuste numero de semanas do ajuste, Pop[cid]
#         vetor indexado por cidade das populacoes, beta e gama parametros do modelo
# retorna: tempo$tempo???(Inf ou Rec)[[cid]][[indice intervalo]]
tempoUniforme <- function(S, I, R, cid, nSemAjuste, Pop, beta, gama)
{
  tempoInfAux1 <- vector()
  tempoRecAux1 <- vector()
  tempoInf <- list()
  tempoRec <- list()
  tempo <- list()
  
  nCid <- length(cid)
  
  for(i in 1:nCid)
  {
    tempoInfCid <- list() # tempos de infeccao
    tempoRecCid <- list() # tempos de recuperacao
    for(j in 1:(nSemAjuste-1))
    {
      # infectados
      n <- I[j+1, i] + (R[j+1, i] - I[j, i] - R[j, i]) # caso em que recuperados j+1 != infectados j
      # gerando os n tempos uniformemente no intervalo
      tempoInfAux1 <- runif(n = n, min = j - 1, max = j)
      # propensidade nos extremos
      h0 <- beta * S[j, i] * I[j, i] / Pop[cid[i]]
      h1 <- beta * S[j+1, i] * I[j+1, i] / Pop[cid[i]]
      # transformacao no processo nao homogeneo de taxa linear
      tempoInfAux1 <- transformaTempo(h0 = h0, h1 = h1, t = tempoInfAux1, i = j - 1)
      tempoInfAux1 <- sort(tempoInfAux1)
      tempoInfCid[[j]] <- tempoInfAux1
      
      # recuperados
      n <- R[j+1, i] - R[j, i]
      # gerando os n tempos uniformemente no intervalo
      tempoRecAux1 <- runif(n = n, min = j - 1, max = j)
      # propensidade nos extremos
      h0 <- gama * I[j, i]
      h1 <- gama * I[j+1, i]
      # transformacao no processo nao homogeneo de taxa linear
      tempoRecAux1 <- transformaTempo(h0 = h0, h1 = h1, t = tempoRecAux1, i = j - 1)
      tempoRecAux1 <- sort(tempoRecAux1)
      tempoRecCid[[j]] <- tempoRecAux1
    }
    # armazenando a lista em um vetor indexado por cidade
    tempoInf[[i]] <- tempoInfCid
    tempoRec[[i]] <- tempoRecCid
  }
  
  tempo <- list(tempoInf = tempoInf, tempoRec = tempoRec)
  return(tempo)
}

#--------------------------------------------------------------------------
# recebe: tempoInfTot[[cid]][[interv]] e tempoRecTot[[cid]][[interv]] listas de instantes de 
#         infeccao e recuperacao indexadas por cidade e intervalo, nCid numero de cidades
#         do ajuste, nSemAjuste numero de semanas de ajuste, Infec[sem, cid] matriz com
#         numero de infectados nas observacoes discretas
# devolve: dupla[[cid]][[interv]]$????(tempo, tipoR ou InfIni)
seqTotalReacoes <- function(tempoInfTot, tempoRecTot, nCid, nSemAjuste, Infec)
{
  seq <- list()
  seq2 <- list()
  
  for(i in 1:nCid)
  {
    for(j in 1:(nSemAjuste - 1))
    {
      seq[[j]] <- sequenciaReacoes(ti = j - 1, tf = j, tempoI = tempoInfTot[[i]][[j]], tempoR = tempoRecTot[[i]][[j]], InfIni = Infec[j, i])
    }
    seq2[[i]] <- seq
  }
  return(seq2)
}

#--------------------------------------------------------------------------
# recebe: ti e tf instantes inicial e final do intervalo, tempoI[n1] e tempoR[n2] vetores
#         com tempos de infeccao e recuperacao do intervalo selecionado, InfIni numero de 
#         infectados no instante inicial do intervalo (viabilidade)
# devolve: dupla$?????(tempo, tipoR ou InfIni)
sequenciaReacoes <- function(ti, tf, tempoI, tempoR, InfIni)
{
  n1 <- length(tempoI)
  n2 <- length(tempoR)
  
  tempo <- vector()
  tipoR <- vector()
  seq <- list()
  
  tempo[1] <- ti
  tipoR[1] <- 0
  
  i <- 1
  j <- 1
  
  # ate o fim de um dos vetores
  while(i <= n1 & j <= n2 & InfIni > 0)
  {
    if(tempoI[i] < tempoR[j])
    {
      # uma infeccao ocorreu
      tempo <- c(tempo, tempoI[i])
      tipoR <- c(tipoR, 1)
      InfIni <- InfIni + 1
      i <- i + 1
    }
    
    else
    {
      # uma recuperacao ocorreu
      tempo <- c(tempo, tempoR[j])
      tipoR <- c(tipoR, 2)
      InfIni <- InfIni - 1
      j <- j + 1
    }
  } # fim do laco em que i < n1 e j < n2
  
  while(i <= n1 & InfIni > 0) # caso hajam infeccoes remanescentes
  {
    tempo <- c(tempo, tempoI[i])
    tipoR <- c(tipoR, 1)
    InfIni <- InfIni + 1
    i <- i + 1
  }
  
  while(j <= n2 & InfIni > 0) # caso hajam recuperacoes remanescentes
  {
    tempo <- c(tempo, tempoR[j])
    tipoR <- c(tipoR, 2)
    InfIni <- InfIni - 1
    j <- j + 1
  }
  
  # o ultimo estado se repete ate o tempo final
  tempo <- c(tempo, tf)
  tipoR <- c(tipoR, 0)
  
  seq <- list(tempo = tempo, tipoR = tipoR, InfIni = InfIni)
}

#-----------------------------------------------------------------------------------------
# mostra a evolucao de S I R descrita pelas duplas ao longo do tempo com a finalidade de 
# impressao
geraMapa <- function(tempo, S, I, R, nCid, nSemAjuste)
{
  traj <- list()
  mapa <- list()
  
  for(cid in 1:nCid)
  {
    for(i in 1:(nSemAjuste - 1))
    {
      Susc <- S[i, cid]
      Infe <- I[i, cid]
      Recu <- R[i, cid]
      traj[[i]] <- geraTrajetoria(tempoInf = tempo$tempoInf[[cid]][[i]], tempoRec = tempo$tempoRec[[cid]][[i]], S = Susc, I = Infe, R = Recu, ti = i - 1, tf = i)
    }
    mapa[[cid]] <- traj
  }
  return(mapa)
}


#---------------------------------------------------------------------------------
# exibe a evolucao do estado do sistema (S, I, R) no decorrer do intervalo [ti, tf]
# traj$???(S, I, R ou tempo)[[indice cidade]][[indice intervalo]]
geraTrajetoria <- function(tempoInf, tempoRec, S, I, R, ti, tf)
{
  n1 <- length(tempoInf)
  n2 <- length(tempoRec)
  
  Susc <- vector()
  Infe <- vector()
  Recu <- vector()
  tempo <- vector()
  tipoR <- vector()
  
  Susc[1] <- S
  Infe[1] <- I
  Recu[1] <- R
  tempo[1] <- ti
  tipoR[1] <- 0
  
  i <- 1
  j <- 1
  
  # ate o fim de um dos vetores
  while(i <= n1 & j <= n2)
  {
    if(tempoInf[i] < tempoRec[j])
    {
      # uma infeccao ocorreu
      Susc <- c(Susc, tail(Susc, n = 1), tail(Susc, n = 1) - 1)
      Infe <- c(Infe, tail(Infe, n = 1), tail(Infe, n = 1) + 1)
      Recu <- c(Recu, tail(Recu, n = 1), tail(Recu, n = 1))
      tempo <- c(tempo, tempoInf[i], tempoInf[i])
      tipoR <- c(tipoR, 0, 1)
      i <- i + 1
    }
    
    else
    {
      # uma recuperacao ocorreu
      Susc <- c(Susc, tail(Susc, n = 1), tail(Susc, n = 1))
      Infe <- c(Infe, tail(Infe, n = 1), tail(Infe, n = 1) - 1)
      Recu <- c(Recu, tail(Recu, n = 1), tail(Recu, n = 1) + 1)
      tempo <- c(tempo, tempoRec[j], tempoRec[j])
      tipoR <- c(tipoR, 0, 2)
      j <- j + 1
    }
  } # fim do laco em que i < n1 e j < n2
  
  while(i <= n1) # caso hajam infeccoes remanescentes
  {
    Susc <- c(Susc, tail(Susc, n = 1), tail(Susc, n = 1) - 1)
    Infe <- c(Infe, tail(Infe, n = 1), tail(Infe, n = 1) + 1)
    Recu <- c(Recu, tail(Recu, n = 1), tail(Recu, n = 1))
    tempo <- c(tempo, tempoInf[i], tempoInf[i])
    tipoR <- c(tipoR, 0, 1)
    i <- i + 1
  }
  
  while(j <= n2) # caso hajam recuperacoes remanescentes
  {
    Susc <- c(Susc, tail(Susc, n = 1), tail(Susc, n = 1))
    Infe <- c(Infe, tail(Infe, n = 1), tail(Infe, n = 1) - 1)
    Recu <- c(Recu, tail(Recu, n = 1), tail(Recu, n = 1) + 1)
    tempo <- c(tempo, tempoRec[j], tempoRec[j])
    tipoR <- c(tipoR, 0, 2)
    j <- j + 1
  }
  
  # o ultimo estado se repete ate o tempo final
  Susc <- c(Susc, tail(Susc, n = 1))
  Infe <- c(Infe, tail(Infe, n = 1))
  Recu <- c(Recu, tail(Recu, n = 1))
  tempo <- c(tempo, tf)
  tipoR <- c(tipoR, 0)
  
  traj <- list(S = Susc, I = Infe, R = Recu, tempo = tempo, tipoR = tipoR)
}

#---------------------------------------------------------------------------------------
# plotar graficos com as trajetorias no intervalo definido no mapa (mapa completo)
plotTraj <- function(mapa, cid, intIni, intFim, discretoI)
{
  tempo <- vector()
  Infec <- vector()
  
  for(i in (intIni+1):(intFim))
  {
    tempo <- c(tempo, mapa[[cid]][[i]]$tempo)
    Infec <- c(Infec, mapa[[cid]][[i]]$I)
  }
  plot(c((intIni):(intFim)), discretoI[(intIni+1):(intFim+1), cid], pch=19, col = 2, xlab = "semanas", ylab = "infectados")
  lines(tempo, Infec, col = "blue")
}


#------------------------------------------------------------------------------------------
# calcula o valor a posteriori
# recebe: priori_a, priori_b valores da priori do modelo, nReacoes numero de reacoes totais
#         do tipo ocorridas durante todo o periodo, integral valor da integral obtida na
#         verossimilhanca
# devolve: valor a posteriori
valorPosteriori <- function(priori_a, priori_b, nReacoes, integral)
{
  valor <- rgamma(n = 1, shape = priori_a + nReacoes, rate = priori_b + integral)
  return(valor)
}

#------------------------------------------------------------------------------------------
# recebe: f_a, f_b valores da funcao nas extremidades a e b, a ponto inicial do intervalo
#         t instante de tempo no intervalo [a, b]
# devolve: (numerico) valor de f_ab(t) onde f_ab eh funcao linearizada nos extremos 
linearizaLambda <- function(f_a, f_b, a, t)
{
  return(f_a + (f_b - f_a)*(t - a))
}

#-------------------------------------------------------------------------------------------
# recebe: dupla[[interv]]$????(tempo, tipoR ou InfIni) dupla que caracteriza trajetoria, 
#         theta (theta[1] = alfa, theta[2] = beta) vetor de parametros, discretoS[sem] e 
#         discretoI[sem] vetores de observacoes discretas do sistema, Pop populacao da 
#         cidade, nSemAjuste numero de semanas do ajuste
# devolve: lambdaP[[interv]]$????(lambda ou h) propensidade linearizada ou real dada
#          a trajetoria
calculaLambdaLinearTotal <- function(dupla, theta, discretoS, discretoI, Pop, nSemAjuste)
{
  lambdaP <- list()
  
  for(j in 1:(nSemAjuste - 1))
  {
    lambdaP[[j]] <- calculaLambdaLinear(tempos = dupla[[j]]$tempo, tipoR = dupla[[j]]$tipoR, theta = theta, Inf_a = discretoI[j], Inf_b = discretoI[j+1], Sus_a = discretoS[j], Sus_b = discretoS[j+1], Pop)
  }
  return(lambdaP)
}  

#------------------------------------------------------------------------------------------
# recebe: tempos[n] vetor de tempos de eventos do intervalo em questao, tipoR[n] vetor com
#         a sequencia de eventos no intervalo estudado ordenado por tempo, theta[reacao]
#         vetor de parametros (theta[1] = alfa, theta[2] = beta), Inf_a e Inf_b valores de
#         infectados nos extremos a e b, Sus_a e Sus_b valores de recuperados nos extremos
#         a e b, Pop valor da populacao na cidade estudada
# devolve: l$????(lambda ou h) propensidades linearizada ou real no intervalo estudado
calculaLambdaLinear <- function(tempos, tipoR, theta, Inf_a, Inf_b, Sus_a, Sus_b, Pop)
{
  lambdaL <- vector()
  h <- vector()
  l <- list()
  
  n <- length(tempos)
  S <- Sus_a
  I <- Inf_a
  indice <- floor(tempos[1])
  
  integral_inf <- 0
  integral_rec <- 0
  
  f_a1 <- theta[1] * Sus_a * Inf_a / Pop
  f_b1 <- theta[1] * Sus_b * Inf_b / Pop
  f_a2 <- theta[2] * Inf_a
  f_b2 <- theta[2] * Inf_b
  
  for(i in 1:(n-1))
  {
    reacao <- tipoR[i]
    dt <- tempos[i+1] - tempos[i]
    integral_inf <- integral_inf + dt * S * I / Pop
    integral_rec <- integral_rec + dt * I
    
    if(reacao == 1)
    {
      lambdaL[i] <- linearizaLambda(f_a = f_a1, f_b = f_b1, a = indice, t = tempos[i])
      h[i] <- theta[1] * S * I / Pop
      S <- S - 1
      I <- I + 1
    }
    else if(reacao == 2)
    {
      lambdaL[i] <- linearizaLambda(f_a = f_a2, f_b = f_b2, a = indice, t = tempos[i])
      h[i] <- theta[2] * I
      I <- I - 1
    }
  }
  lambdaL[1] <- 1
  lambdaL[n] <- 1
  h[1] <- 1
  h[n] <- 1
  
  l <- list(lambda = lambdaL, h = h, integral_inf = integral_inf, integral_rec = integral_rec)
  return(l)
}

# recebe: Lambda$?????(lambda, h, integral_inf, integral_rec) lista contendo a funcao 
#         propensidade linearizada(lambda[n]) - o vetor de funcao propensidade real (h[n]), 
#         as integrais no intervalo de infeccao e recuperacao(sem os parametros), tempo[n] 
#         vetor com os instantes dos eventos no intervalo, theta (theta[1] = beta, theta[2] 
#         = gamma) vetor de parametros do modelo, S_a, S_b suscetiveis no inicio e fim do
#         intervalo, I_a, I_b infectados no inicio e fim do intervalo, Pop populacao da 
#         cidade alvo
# devolve: (valor) derivativo de Radon-Nikodym
derivativoRadonNikodym <- function(lambda, tempo, theta, S_a, S_b, I_a, I_b, Pop)
{
  razaohLambda <- lambda$h / lambda$lambda
  produto <- prod(razaohLambda)
  
  soma_h0 <- 0.5 * ((theta[1]* S_a * I_a / Pop + theta[2] * I_a) + (theta[1]* S_b * I_b / Pop + theta[2] * I_b))
  
  valor <- produto * exp(soma_h0 - (theta[1] * lambda$integral_inf + theta[2] * lambda$integral_rec))
  return(valor)
}

# recebe: duplas lista de duplas que descrevem a trajetoria durante todo o periodo, 
#         indCid indice da cidade estudada, tIni e tFim tempos iniciais e finais considerados
#         no plot, iniI valor de I em tIni
# devolve: lista res$? (t ou I) t tempos, I infectados
ajustaInfPlot <- function(duplas, indCid, tIni, tFim, iniI)
{
  res <- list()
  t <- vector()
  I <- vector()
  I[1] <- iniI         # estou preocupado em plotar apenas infectados
  t[1] <- tIni
  
  for(i in (tIni + 1):(tFim - 1))
  {
    n <- length(duplas[[indCid]][[i]]$tempo)
    for(j in 2:n)
    {
      if(duplas[[indCid]][[i]]$tipoR[j] == 1) #infeccao
      {
        t <- c(t, duplas[[indCid]][[i]]$tempo[j], duplas[[indCid]][[i]]$tempo[j])
        I <- c(I, tail(I, n = 1), tail(I, n = 1) + 1)
      }
      
      else if(duplas[[indCid]][[i]]$tipoR[j] == 2) #recuperacao
      {
        t <- c(t, duplas[[indCid]][[i]]$tempo[j], duplas[[indCid]][[i]]$tempo[j])
        I <- c(I, tail(I, n = 1), tail(I, n = 1) - 1)
      }
      
    }
  }
  res <- list(t = t, I = I)
}

##############################################
#  A EXECUCAO PROPRIAMENTE DITA COMECA AQUI  #
##############################################

#########################################################
# plot de notificacoes para a escolha da semana inicial #
#########################################################
#install.packages("stringr")

library(stringr)

alvo <- 44               #68 09 62 66 43 21 45 82       # cidade a analisar

#z <- c("201252", "201301", "201302", "201303", "201304") 
#z <- c("201301", "201302", "201303", "201304", "201305") 
#z <- c("201302", "201303", "201304", "201305", "201306") 
#z <- c("201303", "201304", "201305", "201306", "201307")
z <- c("201304", "201305", "201306", "201307", "201308") 
#z <- c("201305", "201306", "201307", "201308", "201309") 
#z <- c("201306", "201307", "201308", "201309", "201310") 
#z <- c("201307", "201308", "201309", "201310", "201311") 
#z <- c("201308", "201309", "201310", "201311", "201312") 
#z <- c("201309", "201310", "201311", "201312", "201313") 

color1 <- "blue"
color2 <- "green"

cols <- rep(color1, times = 52)
# selecionando o periodo de ajuste
for(i in 1:length(z))
{
  indice <- match(z[i], names(Not[, alvo]))
  cols[indice] <- color2
}

# salvando grafico
pdf(caminhoAjuste, width = 6, height = 4.5)
notPlot <- barplot(Not_m[,alvo], col = cols, names.arg = FALSE, cex.axis = 1.25, las = 1)
text(notPlot[,1][str_sub(names(Not[,alvo]), start = 6) == "0"], -3.7, srt = 90, adj= 1, xpd = TRUE, labels = rownames(Not)[str_sub(names(Not[,alvo]), start = 6) == "0"] , cex=1.25)
text(notPlot[,1][str_sub(names(Not[,alvo]), start = 6) == "5"], -3.7, srt = 90, adj= 1, xpd = TRUE, labels = rownames(Not)[str_sub(names(Not[,alvo]), start = 6) == "5"] , cex=1.25)
text(notPlot[,1][str_sub(names(Not[,alvo]), start = 5) == "01"], -3.7, srt = 90, adj= 1, xpd = TRUE, labels = rownames(Not)[str_sub(names(Not[,alvo]), start = 5) == "01"] , cex=1.25)
dev.off()

# ESCOLHA DAS CIDADES E DA SEMANA INICIAL

# grupo 1
cid <- c(9, 21, 43, 45, 62, 66, 68, 82)
t0 <- c(201307, 201302, 201304, 201306, 201304, 201304, 201308, 201304)

cid <- c(43, 43)
t0 <- c(201304, 201304)

# grupo 2
#cid <- c(1, 8, 12, 25, 26, 35, 61, 74, 75, 78)
#t0 <- c(201304, 201303, 201304, 201303, 201305, 201252, 201301, 201305, 201304, 201308)

# grupo 3
#cid <- c(03, 07, 13, 15, 29, 33, 44, 57, 58, 80)
#t0 <- c(201305, 201305, 201304, 201252, 201301, 201305, 201304, 201305, 201252, 201301)

ind_t0 <- match(t0, rownames(Not_m))   # indice de t0 na base Not_m
Not_m <- Not_m[, cid]
colnames(Not_m) <- Nomes[cid]

#ficticio
#Not_m <- c(100, 149, 222, 332, 494)

###########################################
# gerar parametros de acordo com a priori #
###########################################
# parametros da distribuicao gamma a priori
a1 <- 1.01
a2 <- 1.01
b1 <- 1
b2 <- 1

beta <- vector()              # vetor com sequencia de betas geradas pelo modelo
gamma <- vector()             # vetor com sequencia de gammas geradas pelo modelo

# priori
beta[1] <- rgamma(n = 1, shape = a1, rate = b1)
gamma[1] <- rgamma(n = 1, shape = a2, rate = b2)

# theta como vetor de parametros
theta <- vector(length = 2)
theta[1] <- beta[1]
theta[2] <- gamma[1]

##################################################
# carregando as observacoes semanais dos estados #
##################################################
# X$?(S, I ou R)[semana do ajuste, cidade]
X <- geraDados(Not_m = Not_m, Pop = Pop, nSemAjuste = nSemAjuste + nSemFut, pSuscPop = pSuscPop, ind_t0 = ind_t0, cid = cid)
#ficticio
X <- list()
X$S <- matrix(nrow = 5, ncol = 2)
X$I <- matrix(nrow = 5, ncol = 2)
X$R <- matrix(nrow = 5, ncol = 2)

X$I[1,] <- 100
X$I[2,] <- 149
X$I[3,] <- 222
X$I[4,] <- 332
X$I[5,] <- 494

X$R[1,] <- 0
X$R[2,] <- 100
X$R[3,] <- 249
X$R[4,] <- 471
X$R[5,] <- 803

X$S <- matrix(nrow = 5, ncol = 2, 1000000) - X$I - X$R 
cid <- c(83, 83)
t0 <- c(0, 0)
Pop <- c(Pop, 1000000)


# contagem do numero de infectados
nInfTotal <- 0
for(i in 1:length(cid))
{
  for(j in 1:(nSemAjuste-1))
  {
    n <- X$I[j+1, i] + (X$R[j+1, i] - X$I[j, i] - X$R[j, i]) # caso em que recuperados j+1 != infectados j
    nInfTotal <- nInfTotal + n
  }
}

# contagem do numero de recuperados
nRecTotal <- sum(X$R[nSemAjuste, ]) - sum(X$R[1, ])

# gerando a trajetoria inicial
# tempo$tempo???(Inf ou Rec)[[indice cidade]][[indice pedaco - subtrajetoria]]
tempo <- tempoUniforme(S = X$S, I = X$I, R = X$R, cid = cid, nSemAjuste = nSemAjuste, Pop = Pop, beta = beta[1], gama = gamma[1])

# gerando as duplas que definem uma trajetoria
tempoTot <- vector()
InfTot <- vector()
duplas <- seqTotalReacoes(tempoInfTot = tempo$tempoInf, tempoRecTot = tempo$tempoRec, nCid = length(cid), nSemAjuste = nSemAjuste, Infec = X$I)

# laco das iteracoes do Metropolis-Hastings
for(j in 1:tamMH)
{
  # gerando as sequencias de duplas
  tempoNovo <- tempoUniforme(S = X$S, I = X$I, R = X$R, cid = cid, nSemAjuste = nSemAjuste, Pop = Pop, beta = beta[1], gama = gamma[1])
  duplasNovo <- seqTotalReacoes(tempoInfTot = tempoNovo$tempoInf, tempoRecTot = tempoNovo$tempoRec, nCid = length(cid), nSemAjuste = nSemAjuste, Infec = X$I)
  
  # variaveis de integrais para o calculo da posteriori  
  integral1 <- 0
  integral2 <- 0
  for(i in 1:length(cid))
  {
    for(k in 1:(nSemAjuste - 1))
    {
      # calculo das funcoes de propensidade linearizada
      lambda <- calculaLambdaLinear(tempos = duplas[[i]][[k]]$tempo, tipoR = duplas[[i]][[k]]$tipoR, theta = theta, Inf_a = X$I[k, i], Inf_b = X$I[k+1, i], Sus_a = X$S[k, i], Sus_b = X$S[k+1, i], Pop = Pop[cid[i]])
      lambdaNovo <- calculaLambdaLinear(tempos = duplasNovo[[i]][[k]]$tempo, tipoR = duplasNovo[[i]][[k]]$tipoR, theta = theta, Inf_a = X$I[k, i], Inf_b = X$I[k+1, i], Sus_a = X$S[k, i], Sus_b = X$S[k+1, i], Pop = Pop[cid[i]])
      # derivativo de Radon-Nikodym
      denominador <- derivativoRadonNikodym(lambda, tempo = duplas[[i]][[k]]$tempo, theta = theta, S_a = X$S[k, i], S_b = X$S[k+1, i], I_a = X$I[k, i], I_b = X$I[k+1, i], Pop = Pop[cid[i]])
      numerador <- derivativoRadonNikodym(lambdaNovo, tempo = duplasNovo[[i]][[k]]$tempo, theta = theta, S_a = X$S[k, i], S_b = X$S[k+1, i], I_a = X$I[k, i], I_b = X$I[k+1, i], Pop = Pop[cid[i]])
      
      print(numerador)
      print(denominador)
      
      # probabilidade de aceitacao e teste MH
      A <- min(1, numerador/denominador)
      u <- runif(n = 1, min = 0, max = 1)
      
      if(u <= A) # caso de troca
      { 
        duplas[[i]][[k]] <- duplasNovo[[i]][[k]]
        lambda <- lambdaNovo
      }
      
      integral1 <- integral1 + lambda$integral_inf
      integral2 <- integral2 + lambda$integral_rec
    } # laco de intervalos
  } # laco de cidades
  
  # calculo das posterioris
  betaNovo <- valorPosteriori(priori_a = a1, priori_b = b1, nReacoes = nInfTotal, integral = integral1)
  gammaNovo <- valorPosteriori(priori_a = a2, priori_b = b2, nReacoes = nRecTotal, integral = integral2)
  
  # atualizacao dos soeficientes
  a1 <- a1 + nInfTotal
  a2 <- a2 + nRecTotal
  b1 <- b1 + integral1
  b2 <- b2 + integral2
  
  # armazena as estimativas
  beta <- c(beta, betaNovo)
  gamma <- c(gamma, gammaNovo)
  theta[1] <- betaNovo
  theta[2] <- gammaNovo
} # laco do numero de iteracoes do MH

# plotando trajetorias
#pdf("/home/nakashima/Dropbox/R/extracaoParametros/gillespie33.pdf", width = 6, height = 4.5)
#plot(prePlot$t, prePlot$I, type = "l", col = "red", xlab = "tempo em semanas", ylab = "infectados", main = Nomes[cid[1]])
#lines(prePlotNovo$t, prePlotNovo$I, col = "springgreen4")
#lines(prePlot2$t, prePlot2$I, col = "purple3")
#points(c(0:4), c(X$I[1, 1], X$I[2, 1], X$I[3, 1], X$I[4, 1], X$I[5, 1]), pch = 21, bg = "red1", col = "red4")
#dev.off()

#prePlot2 <- ajustaInfPlot(duplas = duplas, indCid = 1, tIni = 0, tFim = 5, iniI = X$I[1, 1])
#plot(prePlot2$t, prePlot2$I, type = "l", col = "cyan", xlab = "tempo em semanas", ylab = "infectados", main = Nomes[cid[1]])
#lines(prePlot2$t, prePlot2$I, col = "cyan", xlab = "tempo em semanas", ylab = "infectados", main = Nomes[cid[1]], lwd = 2.5)
#points(c(0:4), c(X$I[1, 1], X$I[2, 1], X$I[3, 1], X$I[4, 1], X$I[5, 1]), pch = 21, bg = "mediumorchid2", col = "mediumorchid4")

# plot da posteriori
parametro <- gamma
corPlot <- "blue"

pdf("/home/nakashima/p80/G3/gamma.pdf", width = 8, height = 5.5)
par(mfrow=c(1,1))
lin <- data.frame(x = c(1:length(parametro)), y = parametro)
linm <- lm(y ~ x, data = lin, subset = 2:4)
plot(c(1:length(parametro)), parametro, type = "l", las = 1, main = expression(paste("Estimativas para ", gamma)), xlab = "iterações", ylab = "", col = corPlot, cex.axis = 1.5, cex.main = 1.5, lwd=1)
par(new=TRUE, oma=c(4,0,0,3))
layout(matrix(1:4,2))
plot(y ~ x, data = lin, xlim=c(71,101), ylim = c(0.8, 0.81), type = "l", las = 1, xlab = "", ylab = "", col = corPlot, cex.axis = 1.5, lwd = 3, yaxt="n")
axis(2,at=c(0.8, 0.805, 0.81), las = 1, cex.axis = 1.5)
dev.off()

#######################
# SIMULACAO GILLESPIE #
#######################

#z <- c("201301", "201302", "201303", "201304", "201305") 
#z <- c("201302", "201303", "201304", "201305", "201306") 
#z <- c("201304", "201305", "201306", "201307", "201308") 
#z <- c("201305", "201306", "201307", "201308", "201309") 
#z <- c("201306", "201307", "201308", "201309", "201310") 
#z <- c("201307", "201308", "201309", "201310", "201311") 
#z <- c("201308", "201309", "201310", "201311", "201312") 
#z <- c("201309", "201310", "201311", "201312", "201313") 

# valor da estimativa dos parametros
betaHat <- tail(beta, n = 1)
gammaHat <- tail(gamma, n = 1)

betaHat <- 1.132
gammaHat <- 1.2

nRep <- 100
cidadeEstudada <- 43
#maximo <- 750

stringCid <- toString(cidadeEstudada)
if(cidadeEstudada < 10)
  stringCid <- paste("0", stringCid, sep = "")

matrizI <- matrix(nrow = nRep, ncol = nSemAjuste)
matrizS <- matrix(nrow = nRep, ncol = nSemAjuste)
matrizDeltaI <- matrix(nrow = nRep, ncol = nSemAjuste)
matrizDeltaS <- matrix(nrow = nRep, ncol = nSemAjuste)
contagemI <- matrix(nrow = nRep, ncol = nSemAjuste)
contagemR <- matrix(nrow = nRep, ncol = nSemAjuste)
contagemI[, 1] <- 0
contagemR[, 1] <- 0
IObs <- vector()
SObs <- vector()

indCid <- match(cidadeEstudada, cid)

# observacoes
IObs[1] <- X$I[1, indCid]
SObs[1] <- Pop[cidadeEstudada] - IObs[1]

for(i in 2:(nSemAjuste+nSemFut))
{
  IObs[i] <- X$I[i, indCid]
  SObs[i] <- SObs[i - 1] - IObs[i]
}

#nome1 <- paste("/home/nakashima/graficos/Gillespie", stringCid, "A.pdf", sep = "")

#pdf(nome1, width = 6, height = 4.5)
#plot(c(0:(nSemAjuste-1)), IObs[1:nSemAjuste], col = "black", pch = 19, ylim = c(0, maximo), xaxt = "n", xlab = "semanas", ylab = "", las = 2, cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.25, main = Nomes[cidadeEstudada])
#axis(1, at=0:4, labels=z[1:5], cex.axis = 1.5)

reacoes1 <- vector()
reacoes2 <- vector()

for(k in 1:nRep)
{
  # contadores de reacoes
  nReacoesI <- 0
  nReacoesR <- 0
  
  SProj <- vector()
  IProj <- vector()
  tempoGil <- vector()
  
  # semana inicial
  SProj[1] <- Pop[cidadeEstudada] - X$I[1, indCid]
  IProj[1] <- X$I[1, indCid]
  
  # observacao da diferenca entre o simulado e o observado (na vida real) em cada semana
  t <- 0
  tempoGil[1] <- t
  for(j in 2:(nSemAjuste))
  {
    while(tail(tempoGil, n = 1) < j - 1)
    {
      # calculo do valor pontual das funcoes de propensidade
      a1 <- betaHat * tail(SProj, n = 1) * tail(IProj, n = 1) / Pop[cidadeEstudada]
      a2 <- gammaHat * tail(IProj, n = 1)
      a0 <- a1 + a2
    
      u <- runif(2, 0, 1)
    
      # verificacao de qual reacao ocorreu
      # infeccao
      if(1.0 * a1 / a0 > u[1])
      {
        SProj <- c(SProj, tail(SProj, n = 1), tail(SProj, n = 1) - 1)
        IProj <- c(IProj, tail(IProj, n = 1), tail(IProj, n = 1) + 1)
        nReacoesI <- nReacoesI + 1
      }
      # recuperacao
      else
      {
        SProj <- c(SProj, tail(SProj, n = 1), tail(SProj, n = 1))
        IProj <- c(IProj, tail(IProj, n = 1), tail(IProj, n = 1) - 1)
        nReacoesR <- nReacoesR + 1
      }
    
      # escolha do tempo com dist exponencial (dependente do estado atual)
      tau <- log(1.0 / u[2])/a0
      t <- t + tau
      tempoGil <- c(tempoGil, t, t)
    }
    contagemI[k, j] <- nReacoesI
    contagemR[k, j] <- nReacoesR
  }
  
  ######################
  #     COMPARACAO     #
  ######################
  SProjS <- vector()
  IProjS <- vector()
  
  # observacoes
  IProjS[1] <- IObs[1]
  SProjS[1] <- SObs[1]
  
  j <- 1
  for(i in 2:nSemAjuste)
  {
    while(tempoGil[j] < i - 1)
    {
      j <- j + 1
    }
    IProjS[i] <- IProj[j - 1]
    SProjS[i] <- SProj[j - 1]
  }
  
  ###############################
  # PLOTANDO O SIMULADO
  #lines(tempoGil, IProj, type = "l", col = k)
  matrizI[k,] <- IProjS
  matrizS[k, ] <- SProjS
  matrizDeltaI[k, ] <- IProjS - IObs[1:nSemAjuste]
  matrizDeltaS[k, ] <- SProjS - SObs[1:nSemAjuste]
  
  reacoes1 <- c(reacoes1, nReacoesI)
  reacoes2 <- c(reacoes2, nReacoesR)
} # fim do laço de contagem de repeticoes

#write.table(contagemI, file = "/home/nakashima/contagemI8.csv", sep = "\t")

for(i in 1:5)
{
  contagemI[, i] <- sort(contagemI[, i])
  contagemR[, i] <- sort(contagemR[, i])
  matrizI[, i] <- sort(matrizI[, i])
}


#points(c(0:(nSemAjuste-1)), IObs[1:nSemAjuste], col = "black", pch = 19)
#dev.off()

colnames(matrizI) <- colnames(matrizI, do.NULL = FALSE, prefix = "semana")
colnames(matrizS) <- colnames(matrizI)
colnames(matrizDeltaI) <- colnames(matrizI)
colnames(matrizDeltaS) <- colnames(matrizI)

rownames(matrizI) <- c(1:nRep)
rownames(matrizS) <- c(1:nRep)
rownames(matrizDeltaI) <- c(1:nRep)
rownames(matrizDeltaS) <- c(1:nRep)
View(matrizDeltaI)
View(matrizDeltaS)

# salvar I
nome2 <- paste("/home/nakashima/IGillespie", stringCid, "A.txt", sep = "")
write.table(matrizDeltaI, nome2, sep="\t")


########################
# ODE ##################
########################
library(deSolve)

#cidadeEstudada <- 43
#indCid <- match(cidadeEstudada, cid)


# padronizado
#N <- Pop[cidadeEstudada]
#s <- X$S[1,indCid]/N
#i <- X$I[1,indCid]/N
#r <- X$R[1,indCid]/N
N <- 1000000
s <- 999900/N
i <- 100/N
r <- 0

estBeta <- 0.8
estGamma <- 1.2

state <- c(s = s, i = i, r = r)
parameters <- c(b = estBeta, g = estGamma)

SIR <- function(t, state, parameters)
{
  with(as.list(c(state, parameters)),{
    ds <- - b * s * i
    di <- b * s * i - g * i
    dr <- g * i
    # return the rate of change
    list(c(ds, di, dr))
  })
}

times <- seq(0, 4, by = 0.01)

# resolvendo o sistema de EDO
out <- ode(y = state, times = times, func = SIR, parms = parameters)
#head(out)
# salvando
#write.table(out, file = "/home/nakashima/Dropbox/R/extracaoParametros/out2.txt", row.names = FALSE, col.names = TRUE)

# voltando ao processo extensivo
IntOut <- out
IntOut[, 2] <- out[, 2] * N
IntOut[, 3] <- out[, 3] * N
IntOut[, 4] <- out[, 4] * N
#head(out)
#plot(c(0:(nSemAjuste-1)), X$I[,1], col = "black", pch = 19)
pdf("/home/nakashima/intensive.pdf", width = 6, height = 4.5)
plot(out[,"time"], out[, "i"], type = "l", lwd = 3, col = "blue", xlab = "tempo", ylab = "i", main = "Dinâmica do processo intensivo")
dev.off()

## dist
# graficos da priori
pdf("/home/nakashima/graficosMM/priori.pdf", width = 6, height = 4.5)
randGamma <- rgamma(n = 10000, 1.01, 1)
hist(randGamma, prob = T, breaks = 100, col = "cyan", main = expression(paste("Priori de ", beta, " e ", gamma)), xlab = "valor", ylab = "densidade")
lines(sort(randGamma), y = dgamma(sort(randGamma), shape = 1.01, rate = 1) , col = "red", lty = 1, lwd = 2)
dev.off()

# grafico da posteriori - beta
pdf("/home/nakashima/posterioriBETAensaio1.pdf", width = 6, height = 4.5)
randGamma <- rgamma(n = 10000, a1, b1)
hist(randGamma, probability = T, breaks = 100, col = "green", main = expression(paste("Posteriori de ", beta)), xlab = "valor", ylab = "", xaxt='n', yaxt='n')
lines(sort(randGamma), y = dgamma(sort(randGamma), shape = a1, rate = b1) , col = "red", lty = 1, lwd = 2.5, xaxt='n')
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.25)
dev.off()

# grafico da posteriori - gamma
pdf("/home/nakashima/posterioriGAMMAensaio1.pdf", width = 6, height = 4.5)
randGamma <- rgamma(n = 10000, a2, b2)
hist(randGamma, probability = T, breaks = 100, col = "green", main = expression(paste("Posteriori de ", gamma)), xlab = "valor", ylab = "", xaxt='n', yaxt='n')
lines(sort(randGamma), y = dgamma(sort(randGamma), shape = a2, rate = b2) , col = "red", lty = 1, lwd = 2.5, xaxt='n')
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.25)
dev.off()

pdf("/home/nakashima/posterioriG2.pdf", width = 6, height = 4.5)
randGamma1 <- rgamma(n = 5000, a1, b1)
h <- hist(randGamma1, breaks = 100)
#h <- hist(randGamma1, freq = FALSE, right = FALSE, breaks = 100, col = "green", main = expression(paste("Posteriori de ", beta)), xlab = "valor", ylab = "densidade")
h$counts <- h$counts / sum(h$counts)
plot(h, freq = TRUE, col = "green", main = expression(paste("Posteriori de ", beta)), xlab = "valor", ylab = "densidade")
densidadeGamma1 <- density(h)
plot(densidadeGamma1, type = "l", lwd = 2.5, col = "blue")
dev.off()

pdf("/home/nakashima/posteriori2.pdf", width = 6, height = 4.5)
hist(rgamma(n = 10000, a2, b2), freq = FALSE, breaks = 100, col = "green", main = expression(paste("Posteriori de ", gamma)), xlab = "valor", ylab = "densidade")
dev.off()

#################################################################################################
########################################################################
#######################
# SIMULACAO GILLESPIE #
#######################

#z <- c("201301", "201302", "201303", "201304", "201305") 
#z <- c("201302", "201303", "201304", "201305", "201306") 
#z <- c("201304", "201305", "201306", "201307", "201308") 
#z <- c("201305", "201306", "201307", "201308", "201309") 
#z <- c("201306", "201307", "201308", "201309", "201310") 
#z <- c("201307", "201308", "201309", "201310", "201311") 
#z <- c("201308", "201309", "201310", "201311", "201312") 
#z <- c("201309", "201310", "201311", "201312", "201313") 

for(indCid in 1:length(cid))
{
  nRep <- 100
  cidadeEstudada <- cid[indCid]
  #t0 <- c(201307, 201302, 201304, 201305, 201305, 201304, 201308, 201305)
  #usar t0 ja definido
  
  z <- seq(from = t0[indCid], to = t0[indCid] + nSemAjuste + nSemFut - 1, by = 1)
  
  if(indCid == 4 || indCid == 9)
  {
    z <- c(201252, 201301, 201302, 201303, 201304, 201305, 201306, 201307, 201308)
  }

  dt <- delta
  colunas <- round(1.0/dt * nSemFut + 1)
  matrizI <- matrix(nrow = nRep, ncol = colunas)
  sortMatrizI <- matrix(nrow = nRep, ncol = colunas)
  resumo <- matrix(nrow = 7, ncol = colunas)

  stringCid <- toString(cidadeEstudada)
  if(cidadeEstudada < 10)
  stringCid <- paste("0", stringCid, sep = "")

  IObs <- vector()
  SObs <- vector()

  #indCid <- match(cidadeEstudada, cid)

  # observacoes
  IObs[1] <- X$I[1, indCid]
  SObs[1] <- Pop[cidadeEstudada] - IObs[1]

  for(i in 2:(nSemAjuste + nSemFut))
  {
    IObs[i] <- X$I[i, indCid]
    SObs[i] <- SObs[i - 1] - IObs[i]
  }

  nome1 <- paste("/home/nakashima/p80/G3/GillespieFut", stringCid, ".pdf", sep = "")

  reacoes1 <- vector()
  reacoes2 <- vector()

  for(k in 1:nRep)
  {
     nReacoesI <- 0
     nReacoesR <- 0
  
    SProj <- vector()
    IProj <- vector()
    tempoGil <- vector()
  
    SProj[1] <- X$S[nSemAjuste, indCid]
    IProj[1] <- X$I[nSemAjuste, indCid]
  
    # observacao da diferenca entre o simulado e o observado (na vida real) em cada semana
    t <- nSemAjuste - 1
    tempoGil[1] <- t
    while(tail(tempoGil, n = 1) < (nSemAjuste + nSemFut - 1))
    {
      #a1 <- betaHat * SProj[i - 1] * IProj[i - 1] / Pop[cidadeEstudada]
      #a2 <- IProj[i - 1]
      a1 <- betaHat * tail(SProj, n = 1) * tail(IProj, n = 1) / Pop[cidadeEstudada]
      a2 <- gammaHat * tail(IProj, n = 1)
      a0 <- a1 + a2
      a3 <- tail(beta, n=1) * tail(SProj, n = 1) * tail(IProj, n = 1) / Pop[cidadeEstudada] + tail(gamma, n=1) * tail(IProj, n = 1)
    
      u <- runif(2, 0, 1)
    
      # verificacao de qual reacao ocorreu
      if(1.0 * a1 / a0 > u[1])
      {
        SProj <- c(SProj, tail(SProj, n = 1), tail(SProj, n = 1) - 1)
        IProj <- c(IProj, tail(IProj, n = 1), tail(IProj, n = 1) + 1)
        #SProj[i] <- SProj[i - 1] - 1
        #IProj[i] <- IProj[i - 1] + 1
        nReacoesI <- nReacoesI + 1
      }
      else
      {
        SProj <- c(SProj, tail(SProj, n = 1), tail(SProj, n = 1))
        IProj <- c(IProj, tail(IProj, n = 1), tail(IProj, n = 1) - 1)
        nReacoesR <- nReacoesR + 1
      }
    
      # escolha do tempo com dist exponencial (dependente do estado atual)
      tau <- log(1.0 / u[2])/a0
      t <- t + tau
      tempoGil <- c(tempoGil, t, t)
    }
  
    ###############################
    # PLOTANDO O SIMULADO
    #lines(tempoGil, IProj, type = "l", col = k)
  
    #matrizI[k, 1] <- IProj[1]
    auxiliar <- IProj[1]
    j <- 1
    for(i in 1:length(tempoGil))
    {
      if(tempoGil[i] > nSemAjuste - 1 + j*dt)
      {
        auxiliar <- c(auxiliar, IProj[i])
        #matrizI[k, j+1] <- IProj[i]
        j <- j + 1
      }
    }
    matrizI[k,] <- auxiliar
  } # fim do laço de contagem de repeticoes
  auxT <- seq((nSemAjuste-1), (nSemAjuste + nSemFut - 1), dt)
  colnames(matrizI) <- auxT

  for(i in 1:colunas)
  {
    sortMatrizI[, i] <- sort(matrizI[, i]) 
  }
  colnames(sortMatrizI) <- auxT

  ##########
  # boxplot
  ##########
  ind1q <- round(0.25 * nRep)
  ind3q <- round(0.75 * nRep)
  meio <- round(0.5 * nRep)

  resumo[1,] <- sortMatrizI[1,]     # menor
  resumo[2,] <- sortMatrizI[nRep,]  # maior
  resumo[3,] <- sortMatrizI[ind1q,] # 1q
  resumo[4,] <- sortMatrizI[ind3q,] # 3q

  dist <- sortMatrizI[ind3q,] - sortMatrizI[ind1q,]


  teste1 <- sortMatrizI[ind1q,] - 1.5 * dist
  teste2 <- sortMatrizI[ind3q,] + 1.5 * dist

  for(j in 1:colunas)
  {
    resumo[5, j] <- max(teste1[j], resumo[1, j]) #limInf
    resumo[6, j] <- min(teste2[j], resumo[2, j]) #limSup
  }

  resumo[7,] <- sortMatrizI[meio,] # mediana

  # plot
  pdf(nome1, width = 6, height = 5)
  maximo <- max(resumo, X$I[,indCid])
  plot(c(0:(nSemAjuste+nSemFut-1)), IObs, col = "black", pch = 19, ylim = c(0, maximo), xaxt = "n", xlab = "semanas", ylab = "", las = 2, cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.25, main = Nomes[cidadeEstudada])
  axis(1, at=0:(nSemAjuste + nSemFut - 1), labels=z, cex.axis = 1.5)
  rect(-100, -100, 4, 2*maximo, col = rgb(131/255, 139/255, 139/255, 0.25), lty=0)
  
  polygon(c(auxT, rev(auxT)), c(resumo[2,], rev(resumo[6,])), col = rgb(153/255, 255/255, 255/255, 1), border = rgb(153/255, 255/255, 255/255, 1))
  polygon(c(auxT,rev(auxT)),c(resumo[4,],rev(resumo[6,])),col= rgb(103/255, 204/255, 255/255, 1), border = rgb(103/255, 204/255, 255/255, 1))
  polygon(c(auxT,rev(auxT)),c(resumo[4,],rev(resumo[7,])),col= rgb(0/255, 128/255, 255/255, 1), border = rgb(0/255, 128/255, 255/255, 1))
  polygon(c(auxT,rev(auxT)),c(resumo[7,],rev(resumo[3,])),col= rgb(0/255, 128/255, 255/255, 1), border = rgb(0/255, 128/255, 255/255, 1))
  polygon(c(auxT,rev(auxT)),c(resumo[3,],rev(resumo[5,])),col= rgb(103/255, 204/255, 255/255, 1), border = rgb(103/255, 204/255, 255/255, 1))
  polygon(c(auxT, rev(auxT)), c(resumo[5,], rev(resumo[1,])), col = rgb(153/255, 255/255, 255/255, 1), border = rgb(153/255, 255/255, 255/255, 1))
  
  lines(auxT, resumo[7,], col = rgb(255/255, 0/255, 0/255, 1), lwd = 2.5)
  
  prePlot2 <- ajustaInfPlot(duplas = duplas, indCid = indCid, tIni = 0, tFim = 5, iniI = X$I[1, indCid])
  lines(prePlot2$t, prePlot2$I, col = rgb(255/255, 64/255, 64/255, 1), lwd = 1.5)
  points(c(0:(nSemAjuste+nSemFut-1)), IObs, col = "black", pch = 19)
  dev.off()
}

##############################################
##################
################################
##################
#############################################
###############

# carregando parametros
betaHat <- 1.493252
gammaHat <- 1
partUnid <- 20
delta <- 1.0 / partUnid


pontual <- matrix(nrow = nSemAjuste, ncol = length(cid))

for(indCid in 1:1)#length(cid))
{
  nRep <- 100
  cidadeEstudada <- cid[indCid]
  #maximo <- 1000
  #t0 <- c(201307, 201302, 201304, 201305, 201305, 201304, 201308, 201305)
  #usar t0 ja definido
  
  z <- seq(from = t0[indCid], to = t0[indCid] + nSemAjuste - 1, by = 1)
  #if(indCid == 6)
  #{
  #  z <- c(201252, 201301, 201302, 201303, 201304)
  #}
  
  dt <- delta
  colunas <- round(1.0/dt * (nSemAjuste - 1) + 1)
  matrizI <- matrix(nrow = nRep, ncol = colunas)
  sortMatrizI <- matrix(nrow = nRep, ncol = colunas)
  resumo <- matrix(nrow = 7, ncol = colunas)
  
  stringCid <- toString(cidadeEstudada)
  if(cidadeEstudada < 10)
    stringCid <- paste("0", stringCid, sep = "")
  
  IObs <- vector()
  SObs <- vector()
  
  #indCid <- match(cidadeEstudada, cid)
  
  # observacoes
  IObs[1] <- X$I[1, indCid]
  SObs[1] <- Pop[cidadeEstudada] - IObs[1]
  
  for(i in 2:(nSemAjuste))
  {
    IObs[i] <- X$I[i, indCid]
    SObs[i] <- SObs[i - 1] - IObs[i]
  }
  
  nome1 <- paste("/home/nakashima/Gillespie2Pre", stringCid, ".pdf", sep = "")
  
  reacoes1 <- vector()
  reacoes2 <- vector()
  
  for(k in 1:nRep)
  {
    nReacoesI <- 0
    nReacoesR <- 0
    
    SProj <- vector()
    IProj <- vector()
    tempoGil <- vector()
    
    SProj[1] <- X$S[1, indCid]
    IProj[1] <- X$I[1, indCid]
    
    # observacao da diferenca entre o simulado e o observado (na vida real) em cada semana
    t <- 0
    tempoGil[1] <- t
    #i <- 2
    while(tail(tempoGil, n = 1) < (nSemAjuste - 1))
    {
      #a1 <- betaHat * SProj[i - 1] * IProj[i - 1] / Pop[cidadeEstudada]
      #a2 <- IProj[i - 1]
      a1 <- betaHat * tail(SProj, n = 1) * tail(IProj, n = 1) / Pop[cidadeEstudada]
      a2 <- gammaHat * tail(IProj, n = 1)
      a0 <- a1 + a2
      a3 <- tail(beta, n=1) * tail(SProj, n = 1) * tail(IProj, n = 1) / Pop[cidadeEstudada] + tail(gamma, n=1) * tail(IProj, n = 1)
      
      u <- runif(2, 0, 1)
      
      if(a1 != 0)
      {
        # verificacao de qual reacao ocorreu
        if(1.0 * a1 / a0 > u[1])
        {
          SProj <- c(SProj, tail(SProj, n = 1), tail(SProj, n = 1) - 1)
          IProj <- c(IProj, tail(IProj, n = 1), tail(IProj, n = 1) + 1)
          #SProj[i] <- SProj[i - 1] - 1
          #IProj[i] <- IProj[i - 1] + 1
          nReacoesI <- nReacoesI + 1
        }
        else
        {
          SProj <- c(SProj, tail(SProj, n = 1), tail(SProj, n = 1))
          IProj <- c(IProj, tail(IProj, n = 1), tail(IProj, n = 1) - 1)
          nReacoesR <- nReacoesR + 1
        }
      
        # escolha do tempo com dist exponencial (dependente do estado atual)
        tau <- log(1.0 / u[2])/a0
        t <- t + tau
        tempoGil <- c(tempoGil, t, t)
      }
      if(a1 == 0)
      {
        SProj <- c(SProj, 0)
        IProj <- c(IProj, 0)
        tempoGil <- c(tempoGil, tail(tempoGil, n=1) + dt)
      }
    }
    
    auxiliar <- IProj[1]
    j <- 1
    for(i in 1:length(tempoGil))
    {
      if(tempoGil[i] > j*dt)
      {
        auxiliar <- c(auxiliar, IProj[i])
        #matrizI[k, j+1] <- IProj[i]
        j <- j + 1
      }
    }
    auxiliar <- auxiliar[1:colunas]
    matrizI[k,] <- auxiliar
  } # fim do laço de contagem de repeticoes
  auxT <- seq(0, (nSemAjuste - 1), dt)
  colnames(matrizI) <- auxT
  
  matrizI[is.na(matrizI)] <- 0
  
  for(i in 1:colunas)
  {
    sortMatrizI[, i] <- sort(matrizI[, i]) 
  }
  colnames(sortMatrizI) <- auxT
  
  ##########
  # boxplot
  ##########
  ind1q <- round(0.25 * nRep)
  ind3q <- round(0.75 * nRep)
  meio <- round(0.5 * nRep)
  
  resumo[1,] <- sortMatrizI[1,]     # menor
  resumo[2,] <- sortMatrizI[nRep,]  # maior
  resumo[3,] <- sortMatrizI[ind1q,] # 1q
  resumo[4,] <- sortMatrizI[ind3q,] # 3q
  
  dist <- sortMatrizI[ind3q,] - sortMatrizI[ind1q,]
  
  
  teste1 <- sortMatrizI[ind1q,] - 1.5 * dist
  teste2 <- sortMatrizI[ind3q,] + 1.5 * dist
  
  for(j in 1:colunas)
  {
    resumo[5, j] <- max(teste1[j], resumo[1, j]) #limInf
    resumo[6, j] <- min(teste2[j], resumo[2, j]) #limSup
  }
  
  resumo[7,] <- sortMatrizI[meio,] # mediana
  
  for(k in 1:(nSemAjuste))
  {
    #pontual[indCid, k + 1] <- resumo[7, k * partUnid + 1]
    pontual[k, indCid] <- resumo[7, (k-1) * partUnid + 1]
  }
  
  
  # plot
  pdf(nome1, width = 6, height = 5)
  maximo <- max(resumo, X$I[1:nSemAjuste,indCid])
  plot(c(0:(nSemAjuste-1)), IObs, col = "black", pch = 19, ylim = c(0, maximo), xaxt = "n", xlab = "semanas", ylab = "", las = 2, cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.25, main = Nomes[cidadeEstudada])
  axis(1, at=0:(nSemAjuste - 1), labels=z, cex.axis = 1.5)

  polygon(c(auxT, rev(auxT)), c(resumo[2,], rev(resumo[6,])), col = rgb(153/255, 255/255, 255/255, 1), border = rgb(153/255, 255/255, 255/255, 1))
  polygon(c(auxT,rev(auxT)),c(resumo[4,],rev(resumo[6,])),col= rgb(103/255, 204/255, 255/255, 1), border = rgb(103/255, 204/255, 255/255, 1))
  polygon(c(auxT,rev(auxT)),c(resumo[4,],rev(resumo[7,])),col= rgb(0/255, 128/255, 255/255, 1), border = rgb(0/255, 128/255, 255/255, 1))
  polygon(c(auxT,rev(auxT)),c(resumo[7,],rev(resumo[3,])),col= rgb(0/255, 128/255, 255/255, 1), border = rgb(0/255, 128/255, 255/255, 1))
  polygon(c(auxT,rev(auxT)),c(resumo[3,],rev(resumo[5,])),col= rgb(103/255, 204/255, 255/255, 1), border = rgb(103/255, 204/255, 255/255, 1))
  polygon(c(auxT, rev(auxT)), c(resumo[5,], rev(resumo[1,])), col = rgb(153/255, 255/255, 255/255, 1), border = rgb(153/255, 255/255, 255/255, 1))
  
  lines(auxT, resumo[7,], col = rgb(255/255, 0/255, 0/255, 1), lwd = 2.5)
  
  #prePlot2 <- ajustaInfPlot(duplas = duplas, indCid = indCid, tIni = 0, tFim = 5, iniI = X$I[1, indCid])
  #lines(prePlot2$t, prePlot2$I, col = rgb(255/255, 64/255, 64/255, 1), lwd = 1.5)
  points(c(0:(nSemAjuste-1)), IObs, col = "black", pch = 19)
  dev.off()
}

perAjusteCid <- X$I[1:nSemAjuste,]

diferenca <- (perAjusteCid - pontual)

difPond <- matrix(nrow = nSemAjuste, ncol = length(cid))

for(i in 1:length(cid))
{
  for(j in 1:nSemAjuste)
  {
    difPond[j,i] <- diferenca[j, i] / perAjusteCid[j, i]
  }
}

difPond <- t(difPond)

qualidade <- vector(length = length(cid))

for(i in 1:length(cid))
{
  qualidade[i] <- sum(difPond[i,])/(nSemAjuste - 1)
}

######################################################
#######################################################
# difusao
#######################################################
#######################################################

betaHat <- 1.177
gammaHat <- 0.8567

tau <- 0.0001
delta <- tau

for(indCid in 1:length(cid))
{
  nRep <- 100
  cidadeEstudada <- cid[indCid]
  N <- Pop[cidadeEstudada]
  #maximo <- 1000
  #t0 <- c(201307, 201302, 201304, 201305, 201305, 201304, 201308, 201305)
  #usar t0 ja definido
  
  z <- seq(from = t0[indCid], to = t0[indCid] + nSemAjuste + nSemFut - 1, by = 1)
  
  dt <- delta
  colunas <- round(1.0/dt * nSemFut + 1)
  matrizI <- matrix(nrow = nRep, ncol = colunas)
  sortMatrizI <- matrix(nrow = nRep, ncol = colunas)
  resumo <- matrix(nrow = 7, ncol = colunas)
  
  stringCid <- toString(cidadeEstudada)
  if(cidadeEstudada < 10)
    stringCid <- paste("0", stringCid, sep = "")
  
  IObs <- vector()
  SObs <- vector()
  
  #indCid <- match(cidadeEstudada, cid)
  
  # observacoes
  IObs[1] <- X$I[1, indCid] / N
  SObs[1] <- 1 * pSuscPop - IObs[1]
  
  for(i in 2:(nSemAjuste + nSemFut))
  {
    IObs[i] <- X$I[i, indCid] / N
    SObs[i] <- (SObs[i - 1] - IObs[i])
  }
  
  nome1 <- paste("/home/nakashima/difusao/0001/G2/difusao", stringCid, ".pdf", sep = "")
  
  reacoes1 <- vector()
  reacoes2 <- vector()
  
  for(k in 1:nRep)
  {
    nReacoesI <- 0
    nReacoesR <- 0
    
    SProj <- vector()
    IProj <- vector()
    tempoDif <- vector()
    
    SProj[1] <- X$S[nSemAjuste, indCid] / N
    IProj[1] <- X$I[nSemAjuste, indCid] / N
    
    # observacao da diferenca entre o simulado e o observado (na vida real) em cada semana
    t <- nSemAjuste - 1
    tempoDif[1] <- t
    #i <- 2
    while(tail(tempoDif, n = 1) < (nSemAjuste + nSemFut - 1))
    {
      # calculo da propensidade pontual
      a1 <- betaHat * tail(SProj, n = 1) * tail(IProj, n = 1)
      a2 <- gammaHat * tail(IProj, n = 1)
      
      z1 <- rnorm(n = 1, mean = 0, sd = 1)
      z2 <- rnorm(n = 1, mean = 0, sd = 1)
      
      sNovo <- tail(SProj, n = 1) - tau * a1 - sqrt(tau * a1 / N) * z1
      iNovo <- tail(IProj, n = 1) + tau * (a1 - a2) + sqrt(tau * a1 / N) * z1 - sqrt(tau * a2 / N) * z2
      
      SProj <- c(SProj, sNovo)
      IProj <- c(IProj, iNovo)

      # escolha do tempo com dist exponencial (dependente do estado atual)
      #tau <- log(1.0 / u[2])/a0
      t <- t + tau
      tempoDif <- c(tempoDif, t)
    }
    
    IProj <- IProj[1:length(matrizI[1,])]
    matrizI[k,] <- IProj
  } # fim do laço de contagem de repeticoes
  auxT <- seq((nSemAjuste-1), (nSemAjuste + nSemFut - 1), dt)
  colnames(matrizI) <- auxT
  
  for(i in 1:colunas)
  {
    sortMatrizI[, i] <- sort(matrizI[, i]) 
  }
  colnames(sortMatrizI) <- auxT
  
  ##########
  # boxplot
  ##########
  ind1q <- round(0.25 * nRep)
  ind3q <- round(0.75 * nRep)
  meio <- round(0.5 * nRep)
  
  resumo[1,] <- sortMatrizI[1,]     # menor
  resumo[2,] <- sortMatrizI[nRep,]  # maior
  resumo[3,] <- sortMatrizI[ind1q,] # 1q
  resumo[4,] <- sortMatrizI[ind3q,] # 3q
  
  dist <- sortMatrizI[ind3q,] - sortMatrizI[ind1q,]
  
  
  teste1 <- sortMatrizI[ind1q,] - 1.5 * dist
  teste2 <- sortMatrizI[ind3q,] + 1.5 * dist
  
  for(j in 1:colunas)
  {
    resumo[5, j] <- max(teste1[j], resumo[1, j]) #limInf
    resumo[6, j] <- min(teste2[j], resumo[2, j]) #limSup
  }
  
  resumo[7,] <- sortMatrizI[meio,] # mediana
  
  # plot
  pdf(nome1, width = 6, height = 5)
  maximo <- max(resumo, X$I[,indCid]/N)
  plot(c(0:(nSemAjuste+nSemFut-1)), IObs, col = "black", pch = 19, ylim = c(0, maximo), xaxt = "n", xlab = "semanas", ylab = "", las = 1, cex.axis = 1, cex.lab = 1, cex.main = 1.25, main = Nomes[cidadeEstudada])
  axis(1, at=0:(nSemAjuste + nSemFut - 1), labels=z, cex.axis = 1.5)

  polygon(c(auxT, rev(auxT)), c(resumo[2,], rev(resumo[6,])), col = rgb(153/255, 255/255, 255/255, 1), border = rgb(153/255, 255/255, 255/255, 1))
  polygon(c(auxT,rev(auxT)),c(resumo[4,],rev(resumo[6,])),col= rgb(103/255, 204/255, 255/255, 1), border = rgb(103/255, 204/255, 255/255, 1))
  polygon(c(auxT,rev(auxT)),c(resumo[4,],rev(resumo[7,])),col= rgb(0/255, 128/255, 255/255, 1), border = rgb(0/255, 128/255, 255/255, 1))
  polygon(c(auxT,rev(auxT)),c(resumo[7,],rev(resumo[3,])),col= rgb(0/255, 128/255, 255/255, 1), border = rgb(0/255, 128/255, 255/255, 1))
  polygon(c(auxT,rev(auxT)),c(resumo[3,],rev(resumo[5,])),col= rgb(103/255, 204/255, 255/255, 1), border = rgb(103/255, 204/255, 255/255, 1))
  polygon(c(auxT, rev(auxT)), c(resumo[5,], rev(resumo[1,])), col = rgb(153/255, 255/255, 255/255, 1), border = rgb(153/255, 255/255, 255/255, 1))
  
  lines(auxT, resumo[7,], col = rgb(255/255, 0/255, 0/255, 1), lwd = 2.5)
  
  #prePlot2 <- ajustaInfPlot(duplas = duplas, indCid = indCid, tIni = 0, tFim = 5, iniI = X$I[1, indCid])
  #lines(prePlot2$t, prePlot2$I/N, col = rgb(255/255, 64/255, 64/255, 1), lwd = 1.5)
  rect(-1, -1, 4, 2*maximo, col = rgb(131/255, 139/255, 139/255, 0.25), lty=0)
  points(c(0:(nSemAjuste+nSemFut-1)), IObs, col = "black", pch = 19)
  dev.off()
}

