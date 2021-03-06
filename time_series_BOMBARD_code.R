rm(list = ls())

setwd(dir = '/Users/Thierry/Desktop')

library(tseries)
library(forecast)
library(TSA)
library(ggfortify)

set.seed(123)

# 1. simuation de la s�rie qui suit un processus stochastique de type ARMA d'ordre 1,1, avec n = 200.

simul = arima.sim(model = list(ar=c(0.4), ma=c(0.54)), n = 200)

# 2. importation des donn�es datas.txt et somme des valeurs de x et de simul
(x = read.table('data.txt', sep = ' '))
simul2 = x[,1] + simul

plot(simul2, main = "ARMA1", xlab = "t", ylab = "X")




              ### OBJECTIF : Ajuster un mod�le SARIMA � nos donn�es x + simul


####### STATIONNARITE

kpss.test(simul2) # On rejette H0, la s�rie n'est pas stationnaire

# Pour mod�liser notre s�rie simul2, nous cherchons � obtenir un bruit blanc faible.
# Ce bruit blanc sera obtenu en retirant les composantes saisonni�res et tendanci�res.



######## COMPOSANTE SAISONNIERE

# Nous pouvons utiliser l'ACF empirique et la PACF empirique pour d�terminer les ordres p et q de notre simulation :
acf(simul2) 
# H0 : les observations entre les lignes bleues sont fortement stationnaires et non corr�l�es
# Nous on cherche � utiliser des observations qui ont un lien entre elle, pour pouvoir estimer le futur avec le pass�
# Ici nous pouvons constater que toutes nos observations sont corr�l�es entre elles, sans aucune convergence vers 0.
# L'allure est sinusoidale ici, il y a donc pr�sence de saisonnalit� dans cette s�rie.

# Puisqu'il y a une composante saisonni�re, on peut la visualiser  l'aide d'un periodogramme

spectrum(simul2)
periodog(simul2) # On remarque 2 peaks : on cherche � connaitre l'indice de la fr�quence maximale
serie = simul2

# indice de la fr�quence pr�dominante :
ind.pers <- which.max(spectrum(serie)$spec)
round(1/(spectrum(serie)$freq[ind.pers]) * frequency(serie))
# la fr�quence max est atteinte � la p�riode r = 10

## On cherche maintenant � valider notre p�riode gr�ce au test de saisonnalit� :
Saison.test(serie, 10) # On rejette H0 l'hypoth�se qui dit que pas de saisonnalit� de p�riode r dans la s�rie
# Il y a donc une pr�sence de saisonnalit� de p�riode 10

# On peut tracer chaque p�riode :
seasonplot(serie, s = 10)

# On cherche maintenant � faire disparaitre la composante saisonni�re avec les diff�rences premi�res :
serie1 = diff(serie, lag = 10) # donc r = 10

# On recheck si notre serie a tjr une saisonnalit� :
acf(serie1) # Il n'y a plus de saisonnalit�. Mais il reste peut �tre une tendance



######## COMPOSANTE TENDANCIERE

kpss.test(serie1) # on recheck la stationnarit�, ici on rejette H0 donc la s�rie n'est pas stationnaire

# Existe-t-il une tendance ? On va utilisr les tests des points mont�es et des discordances :
PtMont.test(serie1) # On va garder H0, l'hypoth�se pas de tendance dans la s�rie.
PtDisc.test(serie1) # On rejette H0 l'hypoth�se pas de tendance dans la s�rie. 
# Il y a donc peut �tre une tendance qu'il faut �liminer. On va quand m�me diff�rencier la s�rie
# pour �tre s�r :

# On diff�rence notre s�rie pour faire disparaitre notre tendance :

new_serie1 = diff(serie1)
PtMont.test(new_serie1) # On garde H0, l'hypoth�se suivante : pas de tendance dans la s�rie
PtDisc.test(new_serie1) # On garde H0, l'hypoth�se suivante : pas de tendance dans la s�rie
# La tendance est nulle

# Nous avons donc diff�renci� une seule fois pour faire disparaitre la tendance.
# ordre d = 1

kpss.test(new_serie1) # Notre s�rie est maintenant stationnaire, on a une p-value > 0.05
# On garde H0.

acfMA(new_serie1)
pacf(new_serie1)
eacf(new_serie1)
# AR, c'est en ligne pour l'Extended Auto Correlation
# MA, c'est en colonne

# Il faut checker le croisement d'un "O" en coin de "X", suivi de d�terminants � droite et en bas
# significatifs.
# On n'a pas de solution en coin mais on va s�lectionner des coefs en coin sans "X"
# on peut voir qu'un modele AR = p = 3 et un mod�le MA = 12 peut faire l'affaire, avec r = 10 et d = 1
# La s�rie devrait �tre estim�e avec un mod�le SARIMA(p,d,q)(0,1,0)[r]


# On cherche � d�terminer si les r�sidus de notre s�rie estim�e forment un bruit blanc :
# Pour rappel un bruit blanc est un ARIMA (0.0.0)
# Bruit blanc fort : processus iid et centr� (E(Xt) = 0).
# Bruit blanc faible : Les variables sont simplement non corr�l�es.
# => l'esp�rance est nulle (s�rie centr�e), la variance est constante, et la covariance entre 2 p�riodes est nulle.
# On peut utiliser le test de LB du tsdiag



##### MODELE 1
(sarima1 = Arima(serie, 
                 order = c(3,1,12),
                 seasonal = list(order=c(0,1,0), 
                                 period=10), 
                 method="ML"))

tsdiag(sarima1)
# 
# On peut garder le mod�le et essayer de le simplifier.
# On a un AIC de 586.56, un AICc de 589.72, un BIC de 638.43.

# Graphique de notre s�rie :
plot(serie, type = "l", ylab = " ")

# Graphique de nos r�sidus estim�s par notre mod�le sarima1 :
plot(serie, type = "l", ylab = " ")+
  lines(sarima1$fitted, type = "l", ylab = " ", col = "red")
# Nous pouvons observer que notre mod�le est plut�t bien ajust�.
# Peut-on le simplifier?


##### SIMPLIFICATION DU MODELE :

# on cherche � d�terminer la significativit� de nos coefficients sous une loi gaussienne
coefs = sarima1$coef
Npar = length(coefs)
ect = sqrt(diag(sarima1$var.coef))[-Npar]
test = coefs/ect
abs(test) <= 1.96

# On peut constater quMe la plupart de nos coefficients sont statistiquement = 0, sauf pour quelques ma.
# On a un ma11 suivi d'un ma12 significatif, ce qui fait qu'on ne peut pas retirer 
# ma12 pour ma10 significatif. On risque d'avoir un mod�le moins robuste, avec un AIC plus �lev�.
# Cependant, on peut se contenter d'un mod�le MA12 :




##### MODELE 2
(sarima2 = Arima(serie, 
                 order = c(0,1,12), 
                 seasonal = list(order=c(0,1,0),
                                 period=10),
                 method="ML"))
tsdiag(sarima2) # La plupart de nos r�sidus pour l'ACF sont = 0, sauf pour le lag = 10.
# Ljung Box a une seule valeur diff�rente de 0, c'est pour le lag = 10.
# On peut garder le mod�le et essayer de le simplifier.
# On a un AIC de 587.61 (plus �lev�), un AICc de 589.69 (plus faible) et un BIC de 629.75 (plus faible).

# Graphique de notre s�rie :
plot(serie, type = "l", ylab = " ")
# Graphique de nos r�sidus estim�s par notre mod�le sarima2 :
plot(serie, type = "l", ylab = " ")+
  lines(sarima2$fitted, type = "l", ylab = " ", col = "red")
# Ce deuxi�me mod�le ajuste tout aussi bien notre s�rie.


# On va voir si nos deux mod�les sont �quivalents :
MI = sarima1 # Mod�le initial
MS = sarima2 # Mod�le simplifi�

nS = length(MS$coef) # nb coefs du mod�le simplifi�
nI = length(MI$coef) # nb coefs du mod�le initial

# On test : H0 : les mod�les sont �quivalents, on s�lectionne le simplifi�
# vs H1 : les mod�les ne sont pas �quivalents, on s�lectionne le mod�le initial
Tobs = -2 * (MS$loglik - MI$loglik)
(pval = 1 - pchisq(Tobs, df = nI - nS))
# Au seuil de 5%, on garde H0, nos deux mod�les sont �quivalents et on selectionne le simplifi�.
# On va se contenter de conserver le sarima2.
# Peut on encore le simplifier ?

coefs2 = sarima2$coef
Npar = length(coefs2)
ect = sqrt(diag(sarima2$var.coef))[-Npar]
test2 = coefs2/ect
abs(test2) <= 1.96
# On ne peut plus simplfier le mod�le, puisque l'on doit retirer ma12, qui est significatif,
# pour baisser le nb de coefficients.



##### MODELE 3 auto.sarima
# On va essayer de comparer les performances de ce mod�le l� avec un modele auto.arima :
(sarima3 = auto.arima(serie, ic = "aic"))
# On a un AIC de 635.01, qui est beaucoup plus �lev� que notre sarima2.
# On pourrait s'interroger sur le fait que notre mod�le n'ai pas introduit de
# composante saisonni�re dans son estimation, alors que l'on avait identifi� des oscillations
# dans l'acf de la serie.
# Du drift a �t� introduit par notre estimation automatique (la moyenne n'est pas nulle).

# Ce mod�le l� a conserv� 3 coefficients, et a utilis� p = 4, d = 1, q = 1.
tsdiag(sarima3) # La plupart de nos r�sidus pour l'ACF sont = 0, sauf pour le lag = 5.
# Ljung Box a une seule valeur diff�rente de 0, c'est pour le lag = 5.
# De plus ce mod�le l� a des r�sidus autocorr�l�s, si on regarde le test d'ind�pendance de LB.



##### PREVISIONS
temps = time(serie)
(previsions = forecast(sarima2, h = 10, conf = 80))

# Dans un intervalle de confiance � 80%, nos pr�visions pour les ann�es 201 et 202 sont
# respectivements �gales � 23.99 et 25.1333, ce qui est conforme � nos bornes inf�rieures et sup�rieures
# de notre intervalle.

# Graphique de nos pr�visions sur une p�riode de 10 :
plot(previsions, 
     plot.conf = T, 
     shaded = T, 
     col = "black", 
     fcol = "blue") 
segments(temps[200],
         serie[200],
         temps[200] + 1/frequency(serie), 
         previsions$mean[1], 
         col = "blue")








