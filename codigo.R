df = read.csv("/home/jamil/Documents/Dump/boxplot-agu.csv")

df = df[2:16,]/df[1:15,] - 1
boxplot(df)
boxplot(df[6:15,])

library(Hmisc)
rcorr(as.matrix(df[6:15,c(1,2)]))
rcorr(as.matrix(df[6:15,c(1,3)]))
rcorr(as.matrix(df[8:15,c(1,4)]))
rcorr(as.matrix(df[6:15,c(2,3)]))
rcorr(as.matrix(df[8:15,c(2,4)]))

df2 = df
df2[11,2] = 0.2641
rcorr(as.matrix(df2[6:15,c(1,2)]))
rcorr(as.matrix(df2[6:15,c(2,3)]))
rcorr(as.matrix(df2[8:15,c(2,4)]))

eigen(cor(df2[6:15,]))

df = read.csv("/home/jamil/Documents/Dump/boxplot-agu.csv")
w = df[16,]

df = log(df[2:16,])-log(df[1:15,])

mu = 0
vAGU = sum(((df$AGU[10:15]-mu)^2)/5)
vPGF = sum(((df$PGF[10:15]-mu)^2)/5)
vPGFN = sum(((df$PGFN[10:15]-mu)^2)/5)
vBC = sum(((df$BACEN[10:15]-mu)^2)/5)

mHighVol = mean(df$PGF[1:5]/2)
vHighVol = var(df$PGF[1:5]/2)

mean(df$PGF[1:5])

mLowVol = mean(df$PGFN[3:7])
vLowVol = var(df$PGFN[3:7])

sqrt(vAGU)
sqrt(vPGF)
sqrt(vBC)

K11 = t(matrix(rep(c(0.92,0.84,0.76,0.68,0.6,0.6,0.6,0.6,0.6,0.6),100000),10,10000))
K12 = t(matrix(rep(c(0.97,0.94,0.91,0.88,0.85,0.85,0.85,0.85,0.85,0.85),100000),10,10000))
K13 = t(matrix(rep(rep(1,10),100000),10,10000))
K14 = t(matrix(rep(c(1.04,1.08,1.12,1.16,1.2,1.2,1.2,1.2,1.2,1.2),100000),10,10000))

wK = runif(10000)
pAGUf = 0.44
cAGUf = 0.176+wK*0.095
wAGUf = 0.433+wK*0.111
pPGFf = 0.99
cPGFf = 0.204+(1-wK)*0.095
wPGFf = 0.261+wK*0.111
pBC   = 0.44
cBC  = (0.00005+runif(10000)*0.00011)*w$BACEN
wBC  = 0.22
pPGFN = 0.571
cPGFN = (0.00005+runif(10000)*0.00011)*w$PGFN
wPGFN = 0.286

rINSS = (1670/w$PGF+0.02*rnorm(10000))
pINSS = rINSS*pPGFf
cINSS = rINSS*cPGFf
wINSS = rINSS*wPGFf
pPGFa = (1-rINSS)*pPGFf
cPGFa = (1-rINSS)*cPGFf
wPGFa = (1-rINSS)*wPGFf

sAGU  = exp(apply(matrix(mu + sqrt(vAGU)*rt(100000,5),10000,10),1,cumsum))
sPGF  = exp(apply(matrix(mu + sqrt(vPGF)*rt(100000,5),10000,10),1,cumsum))
sBC = exp(apply(matrix(mu + sqrt(vBC)*rt(100000,5),10000,10),1,cumsum))
sPGFN = exp(apply(matrix(mu + sqrt(vPGFN)*rt(100000,5),10000,10),1,cumsum))


#Rewrite Capital
valAGUf = sAGU*t(pAGUf + K13*(wAGUf+cAGUf))
valPGFf = sPGF*t(pPGFf + K13*(wPGFf+cPGFf))
valBC   = sBC*t(pBC + K13*(wBC+cBC))
valPGFN = sPGFN*t(pPGFN + K13*(wPGFN+cPGFN))

wv1 = c(sqrt(1/2),sqrt(1/2))
wv2 = c(1/sqrt(3),sqrt(2/3))

############################
### CENÀRIOS FUSÃO TOTAL ###
############################


valAtual = valAGUf + valPGFf + valPGFN + valBC
#Cenário Neutro - União Total
tNew1 = wv1[2]*matrix(0 + sqrt(vPGF)*rt(50000,5),10000,5) + wv1[1]*matrix(mHighVol*wv1[1] + sqrt(vHighVol)*rt(50000,4),10000,5)
tNew2 = matrix(0 + sqrt(vPGF)*rt(50000,5),10000,5)
sNew = exp(apply(cbind(tNew1,tNew2),1,cumsum))
valC1 = sNew*t(pAGUf + pPGFf + pPGFN + pBC + K13*(wAGUf+cAGUf + wPGFf+cPGFf + wBC+cBC + wPGFN+cPGFN))
valOrgC1 = valC1
#
sum(sort(colSums(valOrgC1)) - sort(colSums(valAtual))>0)
a = sort(colSums(valOrgC1)) - sort(colSums(valAtual))
a[2500]
a[5000]
a[7500]
plot(log(sort(colSums(1000000000*valOrgC1))), type = "lines", ylab = "Log dos Gastos em 10 anos")
lines(log(sort(colSums(1000000000*valAtual))), col = "red")

#Cenário Otimista - União Total
tNew1 = wv1[2]*matrix(0 + sqrt(vPGF)*rt(50000,5),10000,5) + wv1[1]*matrix(mHighVol*wv1[1] + sqrt(vHighVol)*rt(50000,4),10000,5)
tNew2 = matrix(0 + sqrt(vPGF)*rt(50000,5),10000,5)
sNew = exp(apply(cbind(tNew1,tNew2),1,cumsum))
valC1 = sNew*t(pAGUf + pPGFf + pPGFN + pBC + K11*(wAGUf+cAGUf + wPGFf+cPGFf + wBC+cBC + wPGFN+cPGFN))
valOrgC1 = valC1
#
sum(sort(colSums(valOrgC1)) - sort(colSums(valAtual))>0)
a = sort(colSums(valOrgC1)) - sort(colSums(valAtual))
a[2500]
a[5000]
a[7500]
plot(log(sort(colSums(1000000000*valOrgC1))), type = "lines", ylab = "Log dos Gastos em 10 anos")
lines(log(sort(colSums(1000000000*valAtual))), col = "red")

#Cenário Eficiente - União Total
tNew1 = wv1[2]*matrix(0 + sqrt(vPGF)*rt(50000,5),10000,5) + wv1[1]*matrix(mHighVol*wv1[1] + sqrt(vHighVol)*rt(50000,4),10000,5)
tNew2 = matrix(0 + sqrt(vPGF)*rt(50000,5),10000,5)
sNew = exp(apply(cbind(tNew1,tNew2),1,cumsum))
valC1 = sNew*t(pAGUf + pPGFf + pPGFN + pBC + K12*(wAGUf+cAGUf + wPGFf+cPGFf + wBC+cBC + wPGFN+cPGFN))
valOrgC1 = valC1
#
sum(sort(colSums(valOrgC1)) - sort(colSums(valAtual))>0)
a = sort(colSums(valOrgC1)) - sort(colSums(valAtual))
a[2500]
a[5000]
a[7500]
plot(log(sort(colSums(1000000000*valOrgC1))), type = "lines", ylab = "Log dos Gastos em 10 anos")
lines(log(sort(colSums(1000000000*valAtual))), col = "red")

#Cenário Otimista com Choques - União Total
tNew1 = wv2[2]*matrix(0 + sqrt(vPGF)*rt(50000,5),10000,5) + wv2[1]*matrix(mHighVol*wv1[1] + sqrt(vHighVol)*rt(50000,4),10000,5)
tNew2 = matrix(0 + sqrt(vPGF)*rt(50000,5),10000,5)
sNew = exp(apply(cbind(tNew1,tNew2),1,cumsum))
valC1 = sNew*t(pAGUf + pPGFf + pPGFN + pBC + K11*(wAGUf+cAGUf + wPGFf+cPGFf + wBC+cBC + wPGFN+cPGFN))
valOrgC1 = valC1
#
sum(sort(colSums(valOrgC1)) - sort(colSums(valAtual))>0)
a = sort(colSums(valOrgC1)) - sort(colSums(valAtual))
a[2500]
a[5000]
a[7500]
plot(log(sort(colSums(1000000000*valOrgC1))), type = "lines", ylab = "Log dos Gastos em 10 anos")
lines(log(sort(colSums(1000000000*valAtual))), col = "red")

#Cenário Eficiente com Choques - União Total
tNew1 = wv2[2]*matrix(0 + sqrt(vPGF)*rt(50000,5),10000,5) + wv2[1]*matrix(mHighVol*wv1[1] + sqrt(vHighVol)*rt(50000,4),10000,5)
tNew2 = matrix(0 + sqrt(vPGF)*rt(50000,5),10000,5)
sNew = exp(apply(cbind(tNew1,tNew2),1,cumsum))
valC1 = sNew*t(pAGUf + pPGFf + pPGFN + pBC + K12*(wAGUf+cAGUf + wPGFf+cPGFf + wBC+cBC + wPGFN+cPGFN))
valOrgC1 = valC1
#
sum(sort(colSums(valOrgC1)) - sort(colSums(valAtual))>0)
a = sort(colSums(valOrgC1)) - sort(colSums(valAtual))
a[2500]
a[5000]
a[7500]
plot(log(sort(colSums(1000000000*valOrgC1))), type = "lines", ylab = "Log dos Gastos em 10 anos")
lines(log(sort(colSums(1000000000*valAtual))), col = "red")


############################
## CENÀRIOS AGU+PGFN TOTAL #
############################

#Cenário Otimista
tNew1 = wv1[2]*matrix(0 + sqrt(vAGU)*rt(50000,5),10000,5) + wv1[1]*matrix(mLowVol*wv1[1] + sqrt(vLowVol)*rt(50000,4),10000,5)
tNew2 = matrix(0 + sqrt(vAGU)*rt(50000,5),10000,5)
sNew = exp(apply(cbind(tNew1,tNew2),1,cumsum))
valNew = sNew*t(pAGUf + pPGFN + K11*(wAGUf+cAGUf+wPGFN+cPGFN))
valOrgC2 = valNew + valPGFf + valBC
#
sum(sort(colSums(valOrgC2)) - sort(colSums(valAtual))>0)
a = sort(colSums(valOrgC2)) - sort(colSums(valAtual))
a[2500]
a[5000]
a[7500]
plot(log(sort(colSums(1000000000*valOrgC2))), type = "lines", ylab = "Log dos Gastos em 10 anos")
lines(log(sort(colSums(1000000000*valAtual))), col = "red")

#Cenário Eficiente
tNew1 = wv1[2]*matrix(0 + sqrt(vAGU)*rt(50000,5),10000,5) + wv1[1]*matrix(mLowVol*wv1[1] + sqrt(vLowVol)*rt(50000,4),10000,5)
tNew2 = matrix(0 + sqrt(vAGU)*rt(50000,5),10000,5)
sNew = exp(apply(cbind(tNew1,tNew2),1,cumsum))
valNew = sNew*t(pAGUf + pPGFN + K12*(wAGUf+cAGUf+wPGFN+cPGFN))
valOrgC2 = valNew + valPGFf + valBC
#
sum(sort(colSums(valOrgC2)) - sort(colSums(valAtual))>0)
a = sort(colSums(valOrgC2)) - sort(colSums(valAtual))
a[2500]
a[5000]
a[7500]
plot(log(sort(colSums(1000000000*valOrgC2))), type = "lines", ylab = "Log dos Gastos em 10 anos")
lines(log(sort(colSums(1000000000*valAtual))), col = "red")


#Cenário Otimista com Choques
tNew1 = wv2[2]*matrix(0 + sqrt(vAGU)*rt(50000,5),10000,5) + wv2[1]*matrix(mLowVol*wv1[1] + sqrt(vLowVol)*rt(50000,4),10000,5)
tNew2 = matrix(0 + sqrt(vAGU)*rt(50000,5),10000,5)
sNew = exp(apply(cbind(tNew1,tNew2),1,cumsum))
valNew = sNew*t(pAGUf + pPGFN + K11*(wAGUf+cAGUf+wPGFN+cPGFN))
valOrgC2 = valNew + valPGFf + valBCx
#
sum(sort(colSums(valOrgC2)) - sort(colSums(valAtual))>0)
a = sort(colSums(valOrgC2)) - sort(colSums(valAtual))
a[2500]
a[5000]
a[7500]
plot(log(sort(colSums(1000000000*valOrgC2))), type = "lines", ylab = "Log dos Gastos em 10 anos")
lines(log(sort(colSums(1000000000*valAtual))), col = "red")

#Cenário Eficiente com Choques
tNew1 = wv2[2]*matrix(0 + sqrt(vAGU)*rt(50000,5),10000,5) + wv2[1]*matrix(mLowVol*wv1[1] + sqrt(vLowVol)*rt(50000,4),10000,5)
tNew2 = matrix(0 + sqrt(vAGU)*rt(50000,5),10000,5)
sNew = exp(apply(cbind(tNew1,tNew2),1,cumsum))
valNew = sNew*t(pAGUf + pPGFN + K12*(wAGUf+cAGUf+wPGFN+cPGFN))
valOrgC2 = valNew + valPGFf + valBC
#
sum(sort(colSums(valOrgC2)) - sort(colSums(valAtual))>0)
a = sort(colSums(valOrgC2)) - sort(colSums(valAtual))
a[2500]
a[5000]
a[7500]
plot(log(sort(colSums(1000000000*valOrgC2))), type = "lines", ylab = "Log dos Gastos em 10 anos")
lines(log(sort(colSums(1000000000*valAtual))), col = "red")

############################
### CENÀRIOS PGF-INSS ######
############################

#Otimista
sINSS = exp(apply(wv1[1]*matrix(0 + sqrt(vBC)*rt(50000,5),10000,10) + wv1[2]*matrix(0 + sqrt(vPGFN)*rt(50000,5),10000,10),1,cumsum))
sPGFa = exp(apply(wv1[1]*matrix(0 + sqrt(vAGU)*rt(50000,5),10000,10) + wv1[2]*matrix(0 + sqrt(vPGF)*rt(50000,5),10000,10),1,cumsum))
valINSS = sINSS*t(pINSS + K13*(wINSS+cINSS))
valPGFa = sPGFa*t(pPGFa + K13*(wPGFa+cPGFa))
valOrgC3 = valAGUf + valPGFa + valPGFN + valBC + valINSS
#
sum(sort(colSums(valOrgC3)) - sort(colSums(valAtual))>0)
a = sort(colSums(valOrgC3)) - sort(colSums(valAtual))
a[2500]
a[5000]
a[7500]
plot(log(sort(colSums(1000000000*valOrgC3))), type = "lines", ylab = "Log dos Gastos em 10 anos")
lines(log(sort(colSums(1000000000*valAtual))), col = "red")

#Deu ruim
sINSS = exp(apply(wv1[1]*matrix(0 + sqrt(vBC)*rt(50000,5),10000,10) + wv1[2]*matrix(0 + sqrt(vPGFN)*rt(50000,5),10000,10),1,cumsum))
sPGFa = exp(apply(wv1[1]*matrix(0 + sqrt(vAGU)*rt(50000,5),10000,10) + wv1[2]*matrix(0 + sqrt(vPGF)*rt(50000,5),10000,10),1,cumsum))
valINSS = sINSS*t(pINSS + K14*(wINSS+cINSS))
valPGFa = sPGFa*t(pPGFa + K14*(wPGFa+cPGFa))
valOrgC3 = valAGUf + valPGFa + valPGFN + valBC + valINSS
#
sum(sort(colSums(valOrgC3)) - sort(colSums(valAtual))>0)
a = sort(colSums(valOrgC3)) - sort(colSums(valAtual))
a[2500]
a[5000]
a[7500]
plot(log(sort(colSums(1000000000*valOrgC3))), type = "lines", ylab = "Log dos Gastos em 10 anos")
lines(log(sort(colSums(1000000000*valAtual))), col = "red")


#End of PGF
muAJ = -0.25
vAJ = 0.005
PGFvar = matrix(muAJ + sqrt(vAJ)*rnorm(100000),10000,10)
sAGUa = exp(apply((matrix(mu + sqrt(vPGF)*rt(100000,5),10000,10)),1,cumsum))
sPGFe = exp(apply(PGFvar,1,cumsum))
sAGUe = t(1+matrix(rep(runif(10000),10),10000,10))*(1-sPGFe)
sAGU = sAGUa + sAGUe

valAGUa = sAGU*t(pAGUf + K12*(wAGUf+cAGUf))
valPGFe = sPGFe*t(pPGFa + K13*(wPGFa+cPGFa))
valINSS = sINSS*t(pINSS + K14*(wINSS+cINSS))
valOrgC4 = valAGUa + valPGFe + valPGFN + valBC + valINSS
#
sum(sort(colSums(valOrgC4)) - sort(colSums(valAtual))>0)
a = sort(colSums(valOrgC4)) - sort(colSums(valAtual))
a[2500]
a[5000]
a[7500]

plot(log(sort(colSums(1000000000*valAtual))), col = "red", type = "lines", ylab = "Log dos Gastos em 10 anos")
lines(log(sort(colSums(1000000000*valOrgC4))), col = "black")
lines(log(sort(colSums(1000000000*valOrgC3))), col = "blue")

valOrgC1 = valC1
valOrgC2 = valNew + valPGFf + valBC
valOrgC3 = valAGUf + valPGFa + valPGFN + valBC + valINSS
valOrgC4 = valAGUe + valPGFe + valPGFN + valBC + valINSS
valAtual = valAGUf + valPGFf + valPGFN + valBC

#valPGFf = sPGF*t(1.31 + K13*(0.30+cPGFf))
#head(tNew1)

sum(sort(colSums(valOrgC1)) - sort(colSums(valAtual))>0)
sum(sort(colSums(valOrgC2)) - sort(colSums(valAtual))>0)
sum(sort(colSums(valOrgC3)) - sort(colSums(valAtual))>0)
sum(sort(colSums(valOrgC4)) - sort(colSums(valAtual))>0)

#PGFN
#K11-metade: 96.85%
#K12-metade: 99.90%

#sort(colSums(valOrgC3))[2500]
#sort(colSums(valAtual))[2500]
#head(t(valOrgC2))
a = sort(colSums(valOrgC4)) - sort(colSums(valAtual))
a[5000]
#IQR

#plot(sort(colSums(log(1000000000*valOrgC1))), type = "lines", ylab = "Log dos Gastos em 10 anos")
#lines(sort(colSums(log(1000000000*valAtual))), col = "red")

#quantile(colSums(1000000000*valOrgC1),0.75)/quantile(colSums(1000000000*valOrgC1),0.25)
#quantile(colSums(1000000000*valOrgC2),0.75)/quantile(colSums(1000000000*valOrgC2),0.25)
#quantile(colSums(1000000000*valAtual),0.75)/quantile(colSums(1000000000*valAtual),0.25)

plot(log(sort(colSums(1000000000*valOrgC3))), type = "lines", ylab = "Log dos Gastos em 10 anos")
lines(log(sort(colSums(1000000000*valAtual))), col = "red")

#plot(sort(colSums(log(1000000000*valOrgC2))), type = "lines", ylab = "Log dos Gastos em 10 anos")
#lines(sort(colSums(log(1000000000*valAtual))), col = "red")

#plot(sort(colSums(sNew)), type = "lines")

quantile(valOrgC1[10,],0.5)
quantile(valOrgC2[10,],0.5)
quantile(valAtual[10,],0.5)


library(tidyverse)
servidores <- readRDS("~/Downloads/agu/servidores.rds")
servidores[[1]]$servidor$id
servidores[[1]]$servidor$situacao$descricao
servidores[[1]]$servidor$pessoa$nome


id = rep(0,8985)
for (i in 1:8985){
  id[i] = servidores[[i]]$servidor$id
  id[i] = servidores[[i]]$fichasCargoEfetivo[[1]]$situacaoServidor
}

remuneracao <- readRDS("~/Downloads/agu/remuneracaolocacao.rds")
counting = remuneracao %>% 
  filter(DESCRICAO_CARGO!="Inválido") %>%
  filter(DESCRICAO_CARGO!="Sem informação") %>%
  filter(ANO==2018) %>% 
  #filter(DESCRICAO_CARGO=="PROCURADOR FEDERAL") %>% 
  group_by(DESCRICAO_CARGO) %>% 
  #group_by(Id_SERVIDOR_PORTAL) %>% 
  summarize(custo = sum(`REMUNERAÇÃO BÁSICA BRUTA (R$)`))

remuneracao[,6] = remuneracao[,6] + remuneracao[,8] + remuneracao[,10] + remuneracao[,12] + remuneracao[,14] + remuneracao[,16]# + remuneracao[,38]

A = remuneracao %>% 
  filter(DESCRICAO_CARGO=="ADVOGADO DA UNIAO") %>%
  group_by(ANO) %>% 
  summarize(custo = sum(`REMUNERAÇÃO BÁSICA BRUTA (R$)`, na.rm = TRUE))

B = remuneracao %>% 
  filter(DESCRICAO_CARGO=="PROCURADOR FEDERAL") %>%
  group_by(ANO) %>% 
  summarize(custo = sum(`REMUNERAÇÃO BÁSICA BRUTA (R$)`, na.rm = TRUE))

C = remuneracao %>% 
  group_by(ANO) %>% 
  summarize(custo = sum(`REMUNERAÇÃO BÁSICA BRUTA (R$)`, na.rm = TRUE))

2.2*0.45

A$custo = 0.758*A$custo
B$custo = 0.764*B$custo
A
B
D = A
D$custo = A$custo + B$custo
D

D = A
D$custo = (A$custo + B$custo)/C$custo
D

A = remuneracao %>% 
  filter(DESCRICAO_CARGO=="ADVOGADO DA UNIAO") %>%
  group_by(Id_SERVIDOR_PORTAL) %>% 
  summarize(custo = sum(`REMUNERAÇÃO BÁSICA BRUTA (R$)`, na.rm = TRUE))

B = remuneracao %>% 
  filter(DESCRICAO_CARGO=="PROCURADOR FEDERAL") %>%
  group_by(Id_SERVIDOR_PORTAL) %>% 
  summarize(custo = sum(`REMUNERAÇÃO BÁSICA BRUTA (R$)`, na.rm = TRUE))

C = remuneracao %>% 
  group_by(Id_SERVIDOR_PORTAL) %>% 
  summarize(custo = sum(`REMUNERAÇÃO BÁSICA BRUTA (R$)`, na.rm = TRUE))

A
B

L1 = l1 %>% 
  filter(Cargo=="Inválido") %>% 
  group_by(Ano) %>% 
  summarize(custo = sum(Remuneracao, na.rm = TRUE))

L2 = l1 %>% 
  filter(Cargo=="Sem informação") %>% 
  group_by(Ano) %>% 
  summarize(custo = sum(Remuneracao, na.rm = TRUE))

L1
L2
L = L1$custo+L2$custo
L

A$custo/C$custo
B$custo/C$custo
C$custo[6] - L[6]

pgfn <- readRDS("~/Downloads/pgfn/remuneracao.rds")
library(tidyverse)
#l1 = aggregate(remuneracao[,30],list(remuneracao$DESCRICAO_CARGO, remuneracao$ANO), sum)
l2 = aggregate(pgfn[,30],list(pgfn$DESCRICAO_CARGO), sum)
sum(l2[46,2])/sum(l2[,2])

rt(1,4)
