rm(list = ls())

setwd(dir = '/Users/Thierry/Desktop')

library(tseries)
library(forecast)
library(TSA)
library(ggfortify)

set.seed(123)

# 1. simuation de la série qui suit un processus stochastique de type ARMA d'ordre 1,1, avec n = 200.

simul = arima.sim(model = list(ar=c(0.4), ma=c(0.54)), n = 200)

# 2. importation des données datas.txt et somme des valeurs de x et de simul
(x = read.table('data.txt', sep = ' '))
simul2 = x[,1] + simul

plot(simul2, main = "ARMA1", xlab = "t", ylab = "X")




              ### OBJECTIF : Ajuster un modèle SARIMA à nos données x + simul


####### STATIONNARITE

kpss.test(simul2) # On rejette H0, la série n'est pas stationnaire

# Pour modéliser notre série simul2, nous cherchons à obtenir un bruit blanc faible.
# Ce bruit blanc sera obtenu en retirant les composantes saisonnières et tendancières.



######## COMPOSANTE SAISONNIERE

# Nous pouvons utiliser l'ACF empirique et la PACF empirique pour déterminer les ordres p et q de notre simulation :
acf(simul2) 
# H0 : les observations entre les lignes bleues sont fortement stationnaires et non corrélées
# Nous on cherche à utiliser des observations qui ont un lien entre elle, pour pouvoir estimer le futur avec le passé
# Ici nous pouvons constater que toutes nos observations sont corrélées entre elles, sans aucune convergence vers 0.
# L'allure est sinusoidale ici, il y a donc présence de saisonnalité dans cette série.

# Puisqu'il y a une composante saisonnière, on peut la visualiser  l'aide d'un periodogramme

spectrum(simul2)
periodog(simul2) # On remarque 2 peaks : on cherche à connaitre l'indice de la fréquence maximale
serie = simul2

# indice de la fréquence prédominante :
ind.pers <- which.max(spectrum(serie)$spec)
round(1/(spectrum(serie)$freq[ind.pers]) * frequency(serie))
# la fréquence max est atteinte à la période r = 10

## On cherche maintenant à valider notre période grâce au test de saisonnalité :
Saison.test(serie, 10) # On rejette H0 l'hypothèse qui dit que pas de saisonnalité de période r dans la série
# Il y a donc une présence de saisonnalité de période 10

# On peut tracer chaque période :
seasonplot(serie, s = 10)

# On cherche maintenant à faire disparaitre la composante saisonnière avec les différences premières :
serie1 = diff(serie, lag = 10) # donc r = 10

# On recheck si notre serie a tjr une saisonnalité :
acf(serie1) # Il n'y a plus de saisonnalité. Mais il reste peut être une tendance



######## COMPOSANTE TENDANCIERE

kpss.test(serie1) # on recheck la stationnarité, ici on rejette H0 donc la série n'est pas stationnaire

# Existe-t-il une tendance ? On va utilisr les tests des points montées et des discordances :
PtMont.test(serie1) # On va garder H0, l'hypothèse pas de tendance dans la série.
PtDisc.test(serie1) # On rejette H0 l'hypothèse pas de tendance dans la série. 
# Il y a donc peut être une tendance qu'il faut éliminer. On va quand même différencier la série
# pour être sûr :

# On différence notre série pour faire disparaitre notre tendance :

new_serie1 = diff(serie1)
PtMont.test(new_serie1) # On garde H0, l'hypothèse suivante : pas de tendance dans la série
PtDisc.test(new_serie1) # On garde H0, l'hypothèse suivante : pas de tendance dans la série
# La tendance est nulle

# Nous avons donc différencié une seule fois pour faire disparaitre la tendance.
# ordre d = 1

kpss.test(new_serie1) # Notre série est maintenant stationnaire, on a une p-value > 0.05
# On garde H0.

acfMA(new_serie1)
pacf(new_serie1)
eacf(new_serie1)
# AR, c'est en ligne pour l'Extended Auto Correlation
# MA, c'est en colonne

# Il faut checker le croisement d'un "O" en coin de "X", suivi de déterminants à droite et en bas
# significatifs.
# On n'a pas de solution en coin mais on va sélectionner des coefs en coin sans "X"
# on peut voir qu'un modele AR = p = 3 et un modèle MA = 12 peut faire l'affaire, avec r = 10 et d = 1
# La série devrait être estimée avec un modèle SARIMA(p,d,q)(0,1,0)[r]


# On cherche à déterminer si les résidus de notre série estimée forment un bruit blanc :
# Pour rappel un bruit blanc est un ARIMA (0.0.0)
# Bruit blanc fort : processus iid et centré (E(Xt) = 0).
# Bruit blanc faible : Les variables sont simplement non corrélées.
# => l'espérance est nulle (série centrée), la variance est constante, et la covariance entre 2 périodes est nulle.
# On peut utiliser le test de LB du tsdiag



##### MODELE 1
(sarima1 = Arima(serie, 
                 order = c(3,1,12),
                 seasonal = list(order=c(0,1,0), 
                                 period=10), 
                 method="ML"))

tsdiag(sarima1)
# 
# On peut garder le modèle et essayer de le simplifier.
# On a un AIC de 586.56, un AICc de 589.72, un BIC de 638.43.

# Graphique de notre série :
plot(serie, type = "l", ylab = " ")

# Graphique de nos résidus estimés par notre modèle sarima1 :
plot(serie, type = "l", ylab = " ")+
  lines(sarima1$fitted, type = "l", ylab = " ", col = "red")
# Nous pouvons observer que notre modèle est plutôt bien ajusté.
# Peut-on le simplifier?


##### SIMPLIFICATION DU MODELE :

# on cherche à déterminer la significativité de nos coefficients sous une loi gaussienne
coefs = sarima1$coef
Npar = length(coefs)
ect = sqrt(diag(sarima1$var.coef))[-Npar]
test = coefs/ect
abs(test) <= 1.96

# On peut constater quMe la plupart de nos coefficients sont statistiquement = 0, sauf pour quelques ma.
# On a un ma11 suivi d'un ma12 significatif, ce qui fait qu'on ne peut pas retirer 
# ma12 pour ma10 significatif. On risque d'avoir un modèle moins robuste, avec un AIC plus élevé.
# Cependant, on peut se contenter d'un modèle MA12 :




##### MODELE 2
(sarima2 = Arima(serie, 
                 order = c(0,1,12), 
                 seasonal = list(order=c(0,1,0),
                                 period=10),
                 method="ML"))
tsdiag(sarima2) # La plupart de nos résidus pour l'ACF sont = 0, sauf pour le lag = 10.
# Ljung Box a une seule valeur différente de 0, c'est pour le lag = 10.
# On peut garder le modèle et essayer de le simplifier.
# On a un AIC de 587.61 (plus élevé), un AICc de 589.69 (plus faible) et un BIC de 629.75 (plus faible).

# Graphique de notre série :
plot(serie, type = "l", ylab = " ")
# Graphique de nos résidus estimés par notre modèle sarima2 :
plot(serie, type = "l", ylab = " ")+
  lines(sarima2$fitted, type = "l", ylab = " ", col = "red")
# Ce deuxième modèle ajuste tout aussi bien notre série.


# On va voir si nos deux modèles sont équivalents :
MI = sarima1 # Modèle initial
MS = sarima2 # Modèle simplifié

nS = length(MS$coef) # nb coefs du modèle simplifié
nI = length(MI$coef) # nb coefs du modèle initial

# On test : H0 : les modèles sont équivalents, on sélectionne le simplifié
# vs H1 : les modèles ne sont pas équivalents, on sélectionne le modèle initial
Tobs = -2 * (MS$loglik - MI$loglik)
(pval = 1 - pchisq(Tobs, df = nI - nS))
# Au seuil de 5%, on garde H0, nos deux modèles sont équivalents et on selectionne le simplifié.
# On va se contenter de conserver le sarima2.
# Peut on encore le simplifier ?

coefs2 = sarima2$coef
Npar = length(coefs2)
ect = sqrt(diag(sarima2$var.coef))[-Npar]
test2 = coefs2/ect
abs(test2) <= 1.96
# On ne peut plus simplfier le modèle, puisque l'on doit retirer ma12, qui est significatif,
# pour baisser le nb de coefficients.



##### MODELE 3 auto.sarima
# On va essayer de comparer les performances de ce modèle là avec un modele auto.arima :
(sarima3 = auto.arima(serie, ic = "aic"))
# On a un AIC de 635.01, qui est beaucoup plus élevé que notre sarima2.
# On pourrait s'interroger sur le fait que notre modèle n'ai pas introduit de
# composante saisonnière dans son estimation, alors que l'on avait identifié des oscillations
# dans l'acf de la serie.
# Du drift a été introduit par notre estimation automatique (la moyenne n'est pas nulle).

# Ce modèle là a conservé 3 coefficients, et a utilisé p = 4, d = 1, q = 1.
tsdiag(sarima3) # La plupart de nos résidus pour l'ACF sont = 0, sauf pour le lag = 5.
# Ljung Box a une seule valeur différente de 0, c'est pour le lag = 5.
# De plus ce modèle là a des résidus autocorrélés, si on regarde le test d'indépendance de LB.



##### PREVISIONS
temps = time(serie)
(previsions = forecast(sarima2, h = 10, conf = 80))

# Dans un intervalle de confiance à 80%, nos prévisions pour les années 201 et 202 sont
# respectivements égales à 23.99 et 25.1333, ce qui est conforme à nos bornes inférieures et supérieures
# de notre intervalle.

# Graphique de nos prévisions sur une période de 10 :
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








