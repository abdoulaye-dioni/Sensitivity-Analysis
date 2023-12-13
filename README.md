# mice_imp_extreme
Analyse de sensibilité avec le module mice dans le cas ou les catégories externes sont manquantes.


```{r,message=FALSE, warning=FALSE}
library(parallel)
library(ggplot2)
library(gridExtra)
library(ggthemes)
library(mice)
```


## Simulation

* Cette fonction simule 4 variables indépendantes `X1 ordinale`, `X2 continue`, `X3 catégorielle nominale` et `X4` binaire

* y est dépendante  et binaire

```{r}
# Définition de la fonction de simulation
simulation_no_hierachical <- function(n, n1_catego, n2_catego, betas) {
  x1 <- factor(sample(1:n1_catego, n, replace = TRUE), ordered = TRUE) #x1 avec n1_catego catégories ordonnées
  x2 <- round(rnorm(n, 20, 4)) # x2 avec une distribution normale
  x3 <- factor(sample(1:n2_catego, n, replace = TRUE), ordered = FALSE) #x3 avec n2_catego catégories non ordonnées
  x4 <- factor(rbinom(n, 1, 0.4)) #  x4 avec une probabilité de succès de 0.4
  X <- model.matrix(~ x1 + x2 + x3 + x4) # design matrix 'eta' pour le modèle
  if (length(betas) != ncol(X)) {
    stop("Produit matriciel impossible")
  }
  eta <- as.vector(X %*% betas)
   proba <- plogis(eta)
  y <- as.factor(rbinom(n, 1, proba)) #variable réponse 'y' binaire basée sur le modèle logistique
  id <- 1:length(y) # numéros d'observation
  data_no_hierachical <- data.frame(id = id, y = y, x1 = x1,  x2 = x2,  x3 = x3,  x4 = x4)
  return(data_no_hierachical)
}
```

 Echantillons simulés
 
```{r}
set.seed(10000) # Générateur du nombre aléatoire

nsim <- 4 # # nombre de simulations 
n <- 2000         # Taille de l'échantillon
n1_catego <- 6     # Nombre de catégories de x1
n2_catego <- 3     # Nombre de catégories de x3
betas <- c(3, -1.5, 2, 1, 2, 2.5, -0.2, 4, 2, 3)
# betas <- c(1, -1., -2, 1, 2, 2, -0.3, 3, -2, 4)
#betas <- c(2, -3., 2, 4, 1, -2, -0.4, 6, -2, 8)  # Coefficients pour le modèle
#betas <- c(5, 3., -2, -4, -1, -2, -0.4, 6, -2, -8)

datas <- vector("list", nsim)
for (s in 1:nsim) {
  datas[[s]] <- simulation_no_hierachical(n = n, n1_catego = n1_catego, n2_catego = n2_catego, betas = betas)
}
```
 





## Données manquantes (focus MNAR)




```{r}
generate_mnar <- function(data, id_obs, ord_var, ord_var.mis, A,probaA, B, probaB, seed) {
  set.seed(seed)
  nA <- round(sum(data[[ord_var]] == A) * probaA)
  idA <- sample(data[[id_obs]][data[[ord_var]] == A], size = nA)
  data[[ord_var.mis]][idA] <- NA
  
  nB <- round(sum(data[[ord_var]] == B) * probaB)
  idB <- sample(data[[id_obs]][data[[ord_var]] == B ], size = nB)
  data[[ord_var.mis]][idB] <- NA
  
  return(data)
}

```



```{r}
probaA = 0.5
probaB = 0.5
A = 1
B = 6
seed = 10000

for(s in 1:nsim){
  datas[[s]]$x1.mis <- datas[[s]]$x1
  datas[[s]] <- generate_mnar(
    data = datas[[s]],
    A =A,
    probaA =probaA, 
    id_obs = "id", 
    B= B,
    probaB = probaB,  
    ord_var = "x1", 
    ord_var.mis = "x1.mis" , 
    seed = seed)
  datas[[s]]$x1.imp <- datas[[s]]$x1.mis
}
```





# Imputation mutliple avec le module mice `contexte MAR`



```{r}

M = 5 # nombre d'imputation

cl <- makeCluster(detectCores() - 2) # nombre de cœurs alloué
clusterExport(cl, c("datas", "M")) #  Exporter 'datas' vers les cluster 
mice.imp <- parLapply(cl, 1:nsim, function(s) {
  mice::mice(
    data= datas[[s]][, c("y", "x1.imp", "x2", "x3", "x4")],  
    m = M,         
    maxit = 10,      
    method = c("logreg", "polr", "pmm", "polyreg", "logreg"),  # Méthodes d'imputation 
    print = FALSE  
  )
})

stopCluster(cl) # Arrêter le cluster parallèle
```



# Algorithme pour l'analyse  de sensibilité



```{r algorithme, fig.align = 'center', out.width = "100%"}
knitr::include_graphics("C:/Users/Dioni Abdoulaye/OneDrive - Université Laval/Bureau/Projet Lynne Moore/Automne 2023/mes captures d'ecrans/algo.png")
```

## `Étape 1, 2 et 3` de l'algorithme



```{r}
REGRESSION <- function(mice.model, daT) {
  n <- ncol(model.matrix(~ y + x2 + x3 + x4, daT)[,-1])
  M <- mice.model$m
  data.imp <- ord.regression <- data.miss <- vector("list", length = M)
  betas <-  matrix(NA, nrow = n, ncol = M)
  zeta <-   matrix(NA, nrow = length(unique(daT$x1)) -1, ncol = M)
  id.miss <- as.numeric(rownames(mice.model$imp$x1.imp))
  eta <- matrix(NA, nrow = length(id.miss), ncol = M)
  
  for (m in 1:M) {
    data.imp[[m]] <- cbind(mice::complete(mice.model, m), daT[, c("id", "x1", "x1.mis")])
    #data.imp[[m]]$x1.imp <- ordered(data.imp[[m]]$x1.imp)
    ord.regression[[m]] <- MASS::polr(x1.imp ~ y + x2 + x3 + x4, method = "probit", Hess = TRUE, data = data.imp[[m]])
    betas[, m] <- as.numeric(coef(ord.regression[[m]]))
    zeta[, m] <- as.numeric(ord.regression[[m]]$zeta)
  }
  
  for (m in 1:M) {
    data.miss[[m]] <- data.imp[[m]][data.imp[[m]]$id %in% id.miss, c("y", "x1.imp", "x2", "x3", "x4")]
    eta[, m] <- model.matrix(~ y + x2 + x3 + x4, data.miss[[m]])[, -1] %*% as.matrix(betas[, m])
  }
  
  dimnames(betas) <- list(names(coef(ord.regression[[1]])), paste0("m", 1:M))
  dimnames(zeta) <- list(paste0("k", 1:(length(unique(daT$x1)) - 1)), paste0("m", 1:M))
  dimnames(eta) <- list(rownames(mice.model$imp$x1.imp), paste0("m", 1:M))
  
  return(list(betas = betas, zeta = zeta, eta = eta, data.imp = data.imp))
}

```


```{r}
cl <- makeCluster(detectCores() - 2) # Nombre de cœurs alloués
clusterExport(cl, c("REGRESSION", "mice.imp", "datas")) # export des objets 
out.regression <- parLapply(cl, 1:nsim, function(s) {
  REGRESSION(
    mice.model = mice.imp[[s]],
    daT = datas[[s]]
  )
})

stopCluster(cl) # Arret du cluster
```




## `Étape 4 et 5` de l'algorithme



```{r}
BUILD <- function(eta, zeta, lambda, seed, delta1 = 1) {
  set.seed(seed) # Générateur du nombre aléatoire
  Theta.lattent <- eta + matrix(rnorm(length(eta), 0, delta1), ncol = ncol(eta)) #  matrice de variables latentes + bruit aléatoire basé sur une distribution normale
  zeta.corrected <- zeta + lambda # Ajout du vecteur lambda pour corriger la matrice zeta
  return(list(Theta.lattent = Theta.lattent, zeta.corrected =zeta.corrected))
}
```



```{r}
create_manyLambda <- function(param, n1_catego, A, B) {
  param_data <- t(tidyr::expand_grid(param, -param)) # matrice de combinaisons de 'param'
  Lambda <- matrix(0, nrow = n1_catego - 1, ncol = ncol(param_data)) # Initialisation
  Lambda[A, ] <- param_data[1, ]
  Lambda[B-1, ] <- param_data[2, ]
  Lambda <- as.data.frame(Lambda)
  colnames(Lambda) <- paste0("lambda", 1:ncol(Lambda))
  return(Lambda)
}


manyLambda <- create_manyLambda(param = seq(0, 4, 1),n1_catego = n1_catego, A = A, B = B)
```


```{r}
manyLambda <- manyLambda[,c(1:2,6,7:8,13)]
manyLambda # differents vecteurs de parametres de sensibilité
```



```{r}

cl <- makeCluster(detectCores() - 2) # nombre de coeurs alloués

clusterExport(cl, c("manyLambda", "BUILD","out.regression","nsim")) # export des objets
out.BUILD <- vector("list", nsim)
out.BUILD <- parLapply(cl, 1:nsim, function(s) {
  lapply(1:ncol(manyLambda), function(l) {
    BUILD(eta = out.regression[[s]]$eta, 
          zeta = out.regression[[s]]$zeta,
          lambda = manyLambda[, l], 
          seed = 10000)
  })
})


stopCluster(cl) # Arret du cluster

```





## `Étape 6` de l'algorithme


```{r}
USE <- function(data.imp,THETA, ZETA ){
  mnar <- matrix(NA, nrow = nrow(THETA), ncol = ncol(THETA)) ## Initialisation 
  data.mnar <- data.imp
  for(i in 1:ncol(THETA)){
    # Remplissage de mnar en fonction de THETA et les intervalles définis par ZETA
    mnar[,i][THETA[,i]<= ZETA[,i][1]] <- 1
    mnar[,i][THETA[,i]> ZETA[,i][1] & THETA[,i]<= ZETA[,i][2]] <- 2
    mnar[,i][THETA[,i]> ZETA[,i][2] & THETA[,i]<= ZETA[,i][3]] <- 3
    mnar[,i][THETA[,i]> ZETA[,i][3] & THETA[,i]<= ZETA[,i][4]] <- 4
    mnar[,i][THETA[,i]> ZETA[,i][4] & THETA[,i]<= ZETA[,i][5]] <- 5
    mnar[,i][THETA[,i]> ZETA[,i][5]] <- 6
  }
   mnar <- as.data.frame( mnar)
  for(i in 1:ncol(THETA)){
    data.mnar[[i]]$x1.mnar <- data.mnar[[i]]$x1.mis
    data.mnar[[i]]$x1.mnar <- as.numeric(data.mnar[[i]]$x1.mnar)
     # Remplacer les valeurs manquantes dans x1.mnar par les valeurs correspondantes dans mnar
    data.mnar[[i]]$x1.mnar[is.na(data.mnar[[i]]$x1.mis)] <- mnar[, i]
    data.mnar[[i]]$x1.mnar <- ordered(data.mnar[[i]]$x1.mnar)
  }
  return(data.mnar)
}
```



```{r}
parallel.use <- function(s,M, manyLambda) {
  
 
  out.USE.chunk <- PROPORTION.chunk <- vector("list", ncol(manyLambda)) #  Initialisation 
  
  # Boucle sur les différentes valeurs de lambda (l)
  for (l in 1: ncol(manyLambda)) {
    # Appel de la fonction USE pour générer des données MNAR avec les valeurs de THETA et ZETA
    out.USE.chunk[[l]] <- PROPORTION.chunk[[l]] <- USE(
      data.imp = out.regression[[s]]$data.imp,
      THETA = out.BUILD[[s]][[l]]$Theta.lattent,
      ZETA = out.BUILD[[s]][[l]]$zeta.corrected)
    
    for (m in 1:M) {
      id_subset <- out.USE.chunk[[l]][[m]]$id %in% datas[[s]]$id[which(is.na(datas[[s]]$x1.mis))]
      
      PROPORTION.chunk[[l]][[m]] <- out.USE.chunk[[l]][[m]][id_subset, c("x1", "x1.imp", "x1.mnar")]
    }
  }
  return(list(out.USE.chunk = out.USE.chunk, PROPORTION.chunk = PROPORTION.chunk))
}
```


```{r}
cl <- makeCluster(detectCores() - 2) # nombre de coeur alloué
clusterExport(cl, c("datas", "USE", "out.BUILD", "parallel.use", "out.regression", "manyLambda", "M"))

use_proportion <- parLapply(cl, 1:nsim, function(s){
   parallel.use(s,M = M, manyLambda = manyLambda)
})
stopCluster(cl)
```






```{r}
# Combiner les résultats de manière hiérarchique
out.USE <- vector("list", nsim)
PROPORTION <- vector("list", nsim)
for (s in 1:nsim) {
  out.USE[[s]] <- PROPORTION[[s]] <- vector("list",  ncol(manyLambda))
  for (l in 1: ncol(manyLambda)) {
    out.USE[[s]][[l]] <- vector("list", M)
    PROPORTION[[s]][[l]] <- vector("list", M)
    for (m in 1:M) {
      out.USE[[s]][[l]][[m]] <- use_proportion[[s]]$out.USE.chunk[[l]][[m]]
      PROPORTION[[s]][[l]][[m]] <- use_proportion[[s]]$PROPORTION.chunk[[l]][[m]]
    }
  }
}
```




# Criteres de selection de lambda pour l'analyse de sensibilité


```{r}
parallel_proportion <- function(s,manyLambda) {
  mnar.table <- mar.table <- full.table <- vector("list", ncol(manyLambda)) #  Initialisation 
  # Création de matrices vides pour stocker les proportions MNAR et MAR
  for (l in 1: ncol(manyLambda)) {
    mnar.table[[l]] <- mar.table[[l]] <- matrix(NA, nrow = M, ncol = length(unique(datas[[1]]$x1)), dimnames = list(paste0("m", 1:M), paste0("X1", 1:length(unique(datas[[1]]$x1)))))
  }
  # Calcul des proportions MNAR et MAR pour chaque ensemble d'imputations m
  for (l in 1: ncol(manyLambda)) {
    for (m in 1:M) {
      mar.table[[l]][m, ] <- round((table(PROPORTION[[s]][[l]][[m]]$x1.imp) / sum(is.na(datas[[s]]$x1.mis))) * 100, 2)
      mnar.table[[l]][m, ] <- round((table(PROPORTION[[s]][[l]][[m]]$x1.mnar) / sum(is.na(datas[[s]]$x1.mis))) * 100, 2)
    }
  }
  
  full.table <- round(t(matrix((table(PROPORTION[[s]][[1]][[1]]$x1)), dimnames = list(paste0("X1", 1:length(unique(datas[[s]]$x1))), "FULL.data")))/sum(is.na(datas[[s]]$x1.mis)) * 100, 2)
  
  for (l in 1:ncol(manyLambda)) {
    mar.table[[l]] <- as.data.frame(t(rbind(full.table, mar.table[[l]])))
    mnar.table[[l]] <- as.data.frame(t(rbind(full.table, mnar.table[[l]])))
    mar.table[[l]]$X1 <- mnar.table[[l]]$X1 <- paste0("X1", 1:length(unique(datas[[1]]$x1)))
    mar.table[[l]] <- mar.table[[l]][, c(ncol(mar.table[[l]]), 1:(ncol(mar.table[[l]]) - 1))]
    mnar.table[[l]] <- mnar.table[[l]][, c(ncol(mnar.table[[l]]), 1:(ncol(mnar.table[[l]]) - 1))]
    row.names(mar.table[[l]]) <- row.names(mnar.table[[l]]) <- NULL
  }
  
  return(list(mar = mar.table, mnar = mnar.table))
}

```




```{r}
cl <- makeCluster(detectCores() -2) # nombre de coeur alloué
clusterExport(cl, c("nsim", "datas", "manyLambda", "M", "PROPORTION", "parallel_proportion"))
suppressWarnings(mar.mnar.result <- parLapply(cl, 1:nsim, function(s){
  parallel_proportion(s , manyLambda = manyLambda)
}))
stopCluster(cl) # arret de cluster
```


```{r}
# Calcul de la moyenne des proportions Full data vs MAR
mar.prop.mean <- data.frame("X1" = mar.mnar.result[[1]]$mar[[1]][, 1], Reduce("+", lapply(mar.mnar.result, function(result) result$mar[[1]][, -1])) / nsim)
mar.prop.mean <- cbind(mar.prop.mean[, 1:2], "all.m" = rowMeans(mar.prop.mean[, -c(1:2)]))

# Calcul des moyennes des proportions Full data vs MNAR pour chacun des lambda
mnar.prop.mean <- lapply(1:ncol(manyLambda), function(l) {
  mnar_df <- data.frame("X1" = mar.mnar.result[[1]]$mnar[[l]][, 1], Reduce("+", lapply(mar.mnar.result, function(result) result$mnar[[l]][, -1])) / nsim)
  mnar_df <- cbind(mnar_df[, 1:2], "all.m" = rowMeans(mnar_df[, -c(1:2)]))
  return(mnar_df)
})

```







Afficher les proportions de MAR pour les catégories de `X1`

```{r}
mar.prop.mean
```

Afficher les résultats de MNAR

```{r}
mnar.prop.mean
```

```{r}
PACKAGE <- function(package){
  for(p in package){
    library(p, character.only = TRUE)
  }
}

```



```{r}

couleur <- c("red", "black",  "orange","green", "purple", "pink", "cyan","magenta", "dimgrey", "blue", "darkred",  "darkgreen", "darkblue", "slategray3","springgreen4","navy", "gold","#FF4D00", "#FF9900","#CCFF00","#3300FF" ,"#00B3FF","turquoise4", "tan","slategray1","sienna3","rosybrown4","peachpuff","mistyrose4","lightyellow","firebrick2")

cl <- makeCluster(detectCores() - 2)
invisible(clusterCall(cl, PACKAGE, c("ggplot2", "ggthemes")))
clusterExport(cl, c("mar.prop.mean", "mnar.prop.mean", "manyLambda", "M", "nsim","couleur"))

comparaison <- parLapply(cl, 1:ncol(manyLambda), function(l) {

# title <- paste(paste0("lambda",l), "=", paste(as.vector(manyLambda[, l]), collapse = ", "))

  resulta <- ggplot() +
    geom_line(data = mar.prop.mean, aes(x = X1, y = FULL.data, group = 1, col = couleur[2])) +
    geom_point(data = mar.prop.mean, aes(x = X1, y = FULL.data, col = couleur[2])) +
    geom_line(data = mar.prop.mean, aes(x = X1, y = all.m, group = 1, col = couleur[4])) +
    geom_point(data = mar.prop.mean, aes(x = X1, y = all.m, col = couleur[4])) +
    geom_line(data = mnar.prop.mean[[l]], aes(x = X1, y = all.m, group = 1, col = couleur[5:(ncol(manyLambda) +4)][l])) +
    geom_point(data = mnar.prop.mean[[l]], aes(x = X1, y = all.m, col = couleur[5:(ncol(manyLambda)+4)][l])) +
    labs(
      y = "Proportion"
    ) +
    scale_color_identity(
      name = "Legend",
      breaks = c(couleur[2], couleur[4], couleur[5:(ncol(manyLambda)+4)]),
      labels = c("simule", "mar", paste0("Mnar",1:ncol(manyLambda))),
      guide = "legend"
    ) +
    theme_clean(base_size = 6) +
   # Supposez que 'l' est l'indice de lambda
ggtitle(bquote(lambda[.(l)] == .(paste(as.vector(manyLambda[, l]), collapse = ", "))))+
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(0,75)

})

stopCluster(cl)
```

Afficher les graphiques en une grille

```{r}
grid.arrange(grobs = comparaison, ncol = 2)
```




# Models d'analyse finale `regression logistique`




## Full data `Échantillons sans valeurs manquantes`

```{r,warning=FALSE, message=FALSE}
cl <- makeCluster(detectCores() -2) # nombre de coeur alloué
clusterExport(cl, c("datas"))

# Modèles GLM
Full.mod <- parLapply(cl, datas, function(data) {
  glm(y ~ x1 + x2 + x3 + x4, family = binomial(link = "logit"), data = data)
})

# Estimés et intervalles de confiance
Full.est <- parLapply(cl, Full.mod, function(model) {
  estime <- model$coefficients
  conf_intervals <- confint(model)
  est_df <- data.frame(
    "estime" = estime,
    "Borne.inf" = conf_intervals[, 1],
    "Borne.sup" = conf_intervals[, 2],
    "Largeur.IC" = conf_intervals[, 2] - conf_intervals[, 1]
  )
  rownames(est_df) <- c("Intercept", paste0("X1", 2:6), "X2", "X32", "X33", "X41")
  return(est_df)
})

stopCluster(cl)
```


```{r}
var_estime <- function(estime_ech, estime_mean){
  return((estime_ech - estime_mean)^2)
}

coverage <- function(value, CI.low, CI.upper){
  ifelse(CI.low <= value & CI.upper >= value,1,0)
}
```


```{r}
full_mcar_mar <- function(true, estimate, nsim){
  estime <- var_by_sim <-  couvert <- Largeur.IC <- data.frame(matrix(NA, nrow = nrow(estimate[[1]]),ncol = nsim))
  for(s in 1:nsim){
    estime[,s] <- estimate[[s]]["estime"]
    Largeur.IC[,s] <- estimate[[s]]["Largeur.IC"]
  var_by_sim[,s] <- var_estime(estime_ech = estimate[[s]][,"estime"], estime_mean = rowMeans(sapply(estimate, function(x) x[,"estime"])) )
  couvert[,s] <- coverage(value = true, CI.low = estimate[[s]]["Borne.inf"], CI.upper = estimate[[s]]["Borne.sup"])
  }

  estime_IC <- data.frame("true" = true,
                          "estime" = rowMeans(estime),
                     "Borne.inf" = rowMeans(estime) - qnorm(0.975)*sqrt(rowSums(var_by_sim)/(nsim -1)),
                     "Borne.sup" = rowMeans(estime) + qnorm(0.975)*sqrt(rowSums(var_by_sim)/(nsim -1)),
                     "sd_estime" = sqrt(rowSums(var_by_sim)/(nsim -1)),
                "Largeur.IC" = rowMeans(Largeur.IC),
                "couvert" = rowMeans(couvert),
                "couvert.inf" = rowMeans(couvert)- qnorm(0.975)*sqrt((rowMeans(couvert)*(1-rowMeans(couvert)))/nsim),
                "couvert.sup" = rowMeans(couvert)+ qnorm(0.975)*sqrt((rowMeans(couvert)*(1-rowMeans(couvert)))/nsim))

   rownames(estime_IC ) <- row.names(estimate[[1]])
  return(estime_IC = estime_IC)
}
```

Moyenne et intervalle de confiance des estimés du total des échantillons simulées

```{r}
estisimule <- full_mcar_mar(true = betas, nsim =nsim,  estimate = Full.est)
estisimule <- round(estisimule,2)
estisimule
```


##  Completes cases analysis `Mécanisme MCAR`


```{r}
cl <- makeCluster(detectCores() -2) # nombre de coeurs alloués
clusterExport(cl, c("datas"))

# Modèles GLM
MCAR.mod <- parLapply(cl, datas, function(data) {
  glm(y ~ x1.mis + x2 + x3 + x4,na.action = na.omit,  family = binomial(link = "logit"), data = data)
})

# Estimés et intervalles de confiance
MCAR.est <- parLapply(cl, MCAR.mod, function(model) {
  estime <- model$coefficients
  conf_intervals <- confint(model)
  est_df <- data.frame(
    "estime" = estime,
    "Borne.inf" = conf_intervals[, 1],
    "Borne.sup" = conf_intervals[, 2],
    "Largeur.IC" = conf_intervals[, 2] - conf_intervals[, 1]
  )
  rownames(est_df) <- c("Intercept", paste0("X1", 2:6), "X2", "X32", "X33", "X41")
  return(est_df)
})

stopCluster(cl) # arret du cluster
```


Moyenne et intervalle de confiance des estimés du total des échantillons sous MCAR

```{r}
estimcar <- full_mcar_mar(true = betas,nsim =nsim, estimate = MCAR.est)
estimcar <- round(estimcar,2)
estimcar
```

## Le module mice `Mécanisme MAR`

```{r}
# Fonction pour effectuer les calculs sur chaque groupe de données
calcul_mar <- function(data) {
  MAR.mod <- list()
  MAR.est <- list()
  
  for (s in 1:nsim) {
    model <- with(data[[s]], glm(y ~ x1.imp + x2 + x3 + x4, family = binomial(link = "logit")))
    
    # Extraire les estimations et les intervalles de confiance
    model_summary <- summary(pool(model), conf.int = TRUE)
    est <- data.frame(
      cbind(
        model_summary[, c("estimate", "2.5 %", "97.5 %")],
        "Largeur.IC" = model_summary[, "97.5 %"] - model_summary[, "2.5 %"]
      )
    )
    
    colnames(est) <- c("estime",  "Borne.inf", "Borne.sup","Largeur.IC")
    rownames(est) <-  c("Intercept", paste0("X1", 2:6), "X2", "X32", "X33", "X41")
    
    MAR.mod[[s]] <- model
    MAR.est[[s]] <- est
  }
  
  return(list(MAR.mod = MAR.mod, MAR.est = MAR.est))
}
```


```{r}
cl <- makeCluster(detectCores()-2) # nombre de coeurs alloués
invisible(clusterCall(cl, PACKAGE, c("mice")))

clusterExport(cl, c("calcul_mar","mice.imp", "nsim"))

# Utiliser parLapply pour effectuer le calcul en parallèle
results <- parLapply(cl, 1:nsim, function(s){
  calcul_mar(data = mice.imp)
})
stopCluster(cl) # Arret des clusters

# Combinez les résultats
tout_MAR.mod <- unlist(lapply(results, function(resulta) resulta$MAR.mod), recursive = FALSE)
tout_MAR.est <- unlist(lapply(results, function(resulta) resulta$MAR.est), recursive = FALSE)
```

```{r}
estimar <- full_mcar_mar(true = betas,nsim = nsim, estimate = tout_MAR.est)
estimar <- round(estimar,2)
estimar
```




## Mécanisme MNAR `Analyse de sensibilité`


```{r}
prara.mnar <- function(s, manyLambda) {
  mnar.list <- MNAR.mod <- MNAR.est <- vector("list",  ncol(manyLambda))
  
  for (l in 1: ncol(manyLambda)) {
    mnar.list[[l]] <- mitml::as.mitml.list(out.USE[[s]][[l]])
    MNAR.mod[[l]] <- with(mnar.list[[l]], glm(y ~ x1.mnar + x2 + x3 + x4, family = binomial(link = "logit")))
    
    MNAR.est[[l]] <- as.data.frame(cbind(
      "estime" = summary(pool(MNAR.mod[[l]]), conf.int = TRUE)[,"estimate"],
      "Borne.inf" = summary(pool(MNAR.mod[[l]]), conf.int = TRUE)[,"2.5 %"],
      "Borne.sup" = summary(pool(MNAR.mod[[l]]), conf.int = TRUE)[,"97.5 %"],
      "Largeur.IC" = summary(pool(MNAR.mod[[l]]), conf.int = TRUE)[,"97.5 %"] - summary(pool(MNAR.mod[[l]]), conf.int = TRUE)[,"2.5 %"]
    ))
    
    MNAR.est[[l]]$term <- NULL
    colnames(MNAR.est[[l]]) <- c("estime", "Borne.inf", "Borne.sup" , "Largeur.IC")
    rownames(MNAR.est[[l]]) <- c("Intercept", paste0("X1", 2:6), "X2", "X32", "X33", "X41")
  }
  
  return(MNAR.est)
}
```


```{r}
cl <- makeCluster(detectCores() - 2)
clusterExport(cl, c("nsim", "manyLambda", "out.USE",  "prara.mnar"))
invisible(clusterCall(cl, PACKAGE, c("mice", "mitml")))

tout.MNAR.est <- parLapply(cl, 1:nsim, function(s){
  prara.mnar(s, manyLambda = manyLambda)
})

stopCluster(cl)

```


```{r}
mnar <- function(true, estimate, l,nsim) {
  estime <- var_by_sim <- couvert <- Largeur.IC <- data.frame(matrix(NA, nrow = nrow(estimate[[1]][[1]]), ncol = nsim))
    coverage.percent <- data.frame(matrix(NA, nrow = nrow(estimate[[1]][[1]]),ncol = 2))
  for (s in 1:nsim) {
    estime[, s] <- estimate[[s]][[l]]["estime"]
    Largeur.IC[, s] = estimate[[s]][[l]]["Largeur.IC"]
   var_by_sim[,s] <- var_estime(estime_ech = estimate[[s]][[l]][, "estime"], estime_mean = rowMeans(sapply(estimate, function(x) x[[l]][,"estime"])))
   
   couvert[,s] <- coverage(value = true, CI.low = estimate[[s]][[l]]["Borne.inf"], CI.upper = estimate[[s]][[l]]["Borne.sup"])
  }
  mnar <- data.frame("true" = true, 
                     "estime" = rowMeans(estime),
                     "Borne.inf" = rowMeans(estime) - qnorm(0.975)*sqrt(rowSums(var_by_sim)/(nsim -1)),
                     "Borne.sup" = rowMeans(estime) + qnorm(0.975)*sqrt(rowSums(var_by_sim)/(nsim -1)),
                     "sd_estime" = sqrt(rowSums(var_by_sim)/(nsim -1)),
                "Largeur.IC" = rowMeans(Largeur.IC),
                "couvert" = rowMeans(couvert),
                "couvert.inf" = rowMeans(couvert)- qnorm(0.975)*sqrt((rowMeans(couvert)*(1-rowMeans(couvert)))/nsim),
                "couvert.sup" = rowMeans(couvert)+ qnorm(0.975)*sqrt((rowMeans(couvert)*(1-rowMeans(couvert)))/nsim))
  rownames(mnar) <- rownames(estimate[[1]][[1]])
  return(mnar = mnar)
}

```


Moyenne et intervalle de confiance des estimés du total des échantillons sous MNAR

```{r}
estimnar <- lapply(1: ncol(manyLambda), function(l) {
  round(mnar(true = betas, estimate = tout.MNAR.est, l = l, nsim = nsim),2)
})
estimnar
```

# Resultats

## Graphiques des estimées et leurs intervalle de confiance

```{r}
estime_graph <- function(data, title, labels, colors) {
  p <- ggplot() +
    #geom_line(data = as.data.frame(data), aes(x = rownames(data), y = true, group = 1, col = labels[1])) +
    geom_point(data = as.data.frame(data), aes(x = rownames(data), y = true, color = labels[1], shape = labels[1]), size = 2) +
    geom_point(data = as.data.frame(data), aes(x = rownames(data), y = estime, color = labels[2], shape = labels[2]), size = 2) +
    geom_errorbar(data = as.data.frame(data), aes(x = rownames(data), y = estime, ymin = Borne.inf, ymax = Borne.sup, color = labels[2]), width = 0.2) +
    labs(title = title, x = "Coefficients", y = "Valeur") +
    scale_color_manual(name = "Legend", values = setNames(colors, labels)) +
    scale_shape_manual(name = "Legend", values = setNames(c(17, 19), labels)) +
    theme_clean(base_size = 7) +
    ylim(range(data)[1], range(data)[2])+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  return(p)
}
```

```{r}
graph1 <- estime_graph(data = estisimule, title = " Beta True  vs Beta Simule ", labels = c( "True", "Simule"),colors = couleur[1:2])
graph2 <- estime_graph(data = estimcar, title = " Beta True  vs Beta Mcar ", labels = c( "True", "Mcar"),colors = couleur[c(1,3)])
graph3 <- estime_graph(data = estimar, title = " Beta True  vs Beta Mar ", labels = c( "True", "Mar"),colors = couleur[c(1,4)])
graph <- list(graph1,graph2,graph3)
```


```{r}
graph.mnar <- lapply(1:ncol(manyLambda), function(l) {
  estime_graph(data = estimnar[[l]], title = paste("Beta True vs",paste( "Beta", paste0("Mnar",l))), labels = c("True", paste0("Mnar",l)), colors = c(couleur[1], couleur[5:(ncol(manyLambda)+4)][l]))
})
```




```{r}
grid.arrange(grobs = c(graph, graph.mnar), ncol = 3)
```


```{r}
estime_biais_function <- function(betas, estisimule, estimcar, estimar, estimnar) {
  
  # Resumé des estimés 
  estimation <- data.frame(
     betas,
    estisimule["estime"],
      estimcar["estime"],
       estimar["estime"],
       sapply(1:length(estimnar),function(l) estimnar[[l]]["estime"])
  )
  
  # Variabilite
  sd_estime <- data.frame(
    estisimule["sd_estime"],
      estimcar["sd_estime"],
       estimar["sd_estime"],
       sapply(1:length(estimnar),function(l) estimnar[[l]]["sd_estime"])
  )

  #  Résumé du biais relatif
  biais_relatif <- as.data.frame(matrix(NA, nrow = nrow(estimation), ncol = ncol(estimation)))
  for (i in 1:ncol(estimation)) {
   biais_relatif[[i]] <- 100 * ((estimation[, i] - estimation[, 1]) / estimation[, 1])
  }
   colnames(estimation) <- colnames(biais_relatif) <-  c("True", "Simule","Mcar","Mar", paste0("Mnar", 1:length(estimnar)))
   colnames(sd_estime) <-  colnames(estimation)[-1]
  rownames(estimation) <- rownames(sd_estime)<- row.names(biais_relatif) <- rownames(estisimule)
  
  return(list(estimation = estimation, sd_estime = sd_estime, biais_relatif = biais_relatif))
}


```


```{r}
EstimationBiais <- estime_biais_function(betas = betas, estisimule = estisimule, estimcar = estimcar, estimar = estimar, estimnar = estimnar)
estimation <- round(EstimationBiais$estimation,2)
sd_estime <-  round(EstimationBiais$sd_estime,2)
biais_relatif <- round(EstimationBiais$biais_relatif,2)
```

## estimés 

```{r}
estimation
```




## Sd  des estimés

```{r}
sd_estime
```



## Biais relatifs

```{r}
biais_relatif
```




```{r}
# Créer un dataframe pour les données de biais
df_biais <- data.frame("variable" =rownames(biais_relatif), biais_relatif)
rownames(df_biais)  <- NULL

# Utiliser la fonction gather pour regrouper les données
df_biais_long <- tidyr::gather(df_biais, Method, Bias, -variable)

labels_methods <- unique(df_biais_long$Method)

# Créer le graphique ggplot avec des lignes connectant les points de même couleur
gg <- ggplot(df_biais_long, aes(x = variable, y = Bias, color = Method, group = Method)) +
  geom_point() +
  geom_line() +
  scale_color_manual(
    values = couleur[1:ncol(df_biais)],
    labels = labels_methods,
    breaks = labels_methods
  ) +
  labs(
    x = "Coefficients",
    y = "Relative bias %",
    title = "Comparaison des biais",
    color = "Méthode"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(range(df_biais[-1])[1], range(df_biais[-1])[2])

```



```{r}
# Afficher le graphique
print(gg)
```



