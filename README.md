# mice_imp_extreme
Sensitivity analysis approach in the presence of a  missing not at random ordinal independent variable: non hierachical context


```{r,message=FALSE, warning=FALSE}
library(parallel) # Load the library for executing operations in parallel
library(ggplot2) #  Load the library for creating elegant and informative plots
library(gridExtra) # Load the library for organizing multiple plots in a single layout
library(ggthemes) # Load the library for applying aesthetic themes to ggplot2 plots
library(mice) # Load the library for performing missing data imputation
```


## Simulation

The function "simulation_no_hierachical" simulates a dataset with binary response and predictor variables based on specified characteristics, providing a useful tool for testing and experimenting with statistical models

```{r}
simulation_no_hierachical <- function(n, n1_catego, n2_catego, betas) {
  x1 <- factor(sample(1:n1_catego, n, replace = TRUE), ordered = TRUE) # Generate  x1 with n1_catego ordered categories
  x2 <- round(rnorm(n, 20, 4))  # Generate  x2 with a normal distribution
  x3 <- factor(sample(1:n2_catego, n, replace = TRUE), ordered = FALSE) # Generate  x3 with n2_catego unordered categories
  x4 <- factor(rbinom(n, 1, 0.4)) # Generate x4 with a success probability of 0.4
  X <- model.matrix(~ x1 + x2 + x3 + x4)  # Create the design matrix 'X' for the model
  if (length(betas) != ncol(X)) { 
    stop("Matrix multiplication not possible. Number of betas does not match the number of columns in the design matrix.")
  }
  eta <- as.vector(X %*% betas) # Calculate the linear predictor 'eta' 
  proba <- plogis(eta) # Calculate the probability of success 
  y <- as.factor(rbinom(n, 1, proba)) # Generate binary response variable 'y' 
  id <- 1:length(y)
  data_no_hierarchical <- data.frame(id = id, y = y, x1 = x1, x2 = x2, x3 = x3, x4 = x4) # Create a data frame
  return(data_no_hierarchical)
}
```

 Input Parameters for Data generation
 
```{r}
set.seed(10000) # Random Number Generator.
nsim <- 1000 # # number simulations 
n <- 2000        # The total number of observations in the simulated dataset.
n1_catego <- 6     # The number of categories for the ordered factor variable x1.
n2_catego <- 3     # The number of categories for the unordered factor variable x3.
betas <- c(3, -1.5, 2, 1, 2, 2.5, -0.2, 4, 2, 3) # A vector of coefficients to be used in the linear predictor.
datas <- lapply(1:nsim, function(simulation) {
   simulation_no_hierachical(n = n, n1_catego = n1_catego, n2_catego = n2_catego, betas = betas)
  })
```
 


## Missing data simulation
This function is useful for simulating data with missing values that are not at random.




```{r}
# Function to generate MNAR data
generate_mnar <- function(data, id_obs, ord_var, ord_var.mis, A, probaA, B, probaB, seed) {
  set.seed(seed) # Set the seed for reproducibility
  nA <- round(sum(data[[ord_var]] == A) * probaA) # Calculate the number of missing values for category A based on the specified probability
  idA <- sample(data[[id_obs]][data[[ord_var]] == A], size = nA) #   Sample random IDs from category A to introduce missing values
  data[[ord_var.mis]][idA] <- NA   # Set missing values in the specified variable for category A
  nB <- round(sum(data[[ord_var]] == B) * probaB)  # Calculate the number of missing values for category B based on the specified probability
  idB <- sample(data[[id_obs]][data[[ord_var]] == B ], size = nB) # Sample random IDs from category B to introduce missing values
  data[[ord_var.mis]][idB] <- NA  # Set missing values in the specified variable for category B
  return(data)
}


```


Input Parameters for Missingness generation
```{r}
probaA = 0.5
probaB = 0.5
A = 1
B = 6
seed = 10000

datas <- lapply(datas, function(data) {
  data$x1.mis <- data$x1
  data <- generate_mnar(data = data, A = A, probaA = probaA, id_obs = "id", B = B, probaB = probaB, ord_var = "x1", ord_var.mis = "x1.mis", seed = seed)
  data$x1.imp <- data$x1.mis
  return(data)
})
```





# First step : Perform multiple imputations using the mice package

```{r}
M = 10 # Number of imputations
cl <- makeCluster(round(detectCores()/2)) # Create a parallel cluster with half the available cores
clusterExport(cl, c("datas", "M")) # Export necessary variables to the cluster
mice.imp <- parLapply(cl, 1:nsim, function(s) { # Perform multiple imputations in parallel using parLapply
  mice::mice(
    data = datas[[s]][, c("y", "x1.imp", "x2", "x3", "x4")],  
    m = M,         
    maxit = 10,      
    method = c("logreg", "polr", "pmm", "polyreg", "logreg"),  # Imputation methods 
    print = FALSE  
  )
})
stopCluster(cl) # Stop the parallel cluster
```



# Algorithme


## `Étape 1, 2 et 3` de l'algorithme

This function REGRESSION  perform regression analyses on imputed datasets obtained from the mice package

```{r}
REGRESSION <- function(mice.model, daT) {
  n <- ncol(model.matrix(~ y + x2 + x3 + x4, daT)[,-1])
  M <- mice.model$m
  id.miss <- daT$id[which(is.na(daT$x1.mis))]
  data.imp <- ord.regression <- data.miss <- vector("list", length = M)
  betas <- matrix(NA, nrow = n, ncol = M)
  zeta <- matrix(NA, nrow = length(unique(daT$x1)) - 1, ncol = M)
  eta <- matrix(NA, nrow = length(id.miss), ncol = M)
 data.imp <- lapply(1:M, function(m) { # Generate completed datasets
    completed_data <- mice::complete(mice.model, m)
    selected_columns <- daT[, c("id", "x1", "x1.mis")]
    combined_data <- cbind(completed_data, selected_columns)
    return(combined_data)
  }) 
  
  # Perform ordered logistic regression on each imputed dataset
  ord.regression <- lapply(data.imp, function(imp_data) {
    return(MASS::polr(x1.imp ~ y + x2 + x3  + x4, method = "probit", Hess = TRUE, data = imp_data))
  })
  
  # Extract regression coefficients and zeta values
  betas <- sapply(ord.regression, function(ord_reg) {
    return(as.numeric(coef(ord_reg)))
  })
  
  zeta <- sapply(ord.regression, function(ord_reg) {
    return(as.numeric(ord_reg$zeta))
  })
  
  # Extract relevant columns for missing data points
  data.miss <- lapply(data.imp, function(imp_data) {
    return(imp_data[imp_data$id %in% id.miss, c("y", "x1.imp", "x2", "x3", "x4")])
  })
  
  # Calculate eta values for missing data points
  eta <- sapply(1:M, function(m) {
    model_matrix <- model.matrix(~ y + x2 + x3 + x4, data.miss[[m]])[, -1]
    return(model_matrix %*% betas[, m])
  })
  dimnames(betas) <- list(names(coef(ord.regression[[1]])), paste0("m", 1:M))
  dimnames(zeta) <- list(paste0("k", 1:(length(unique(daT$x1)) - 1)), paste0("m", 1:M))
  dimnames(eta) <- list(as.character(id.miss), paste0("m", 1:M))
  return(list(betas = betas, zeta = zeta, eta = eta, data.imp = data.imp))
}
```


```{r}
cl <- makeCluster(round(detectCores()/2)) # Nombre de cœurs alloués
clusterExport(cl, c("REGRESSION", "mice.imp", "datas")) # export des objets 
out.regression <- parLapply(cl, 1:nsim, function(s) {
  REGRESSION(
    mice.model = mice.imp[[s]],
    daT = datas[[s]]
  )
})
stopCluster(cl) 
```




## `Étape 4 et 5` de l'algorithme

This function named BUILD that appears to generate latent variables (Theta.lattent) and corrects a matrix (zeta) by adding a vector (lambda).

```{r}
BUILD <- function(eta, zeta, lambda, seed, delta1 = 1) {
  set.seed(seed)
 # Generate latent variables by adding random noise to eta
  Theta.lattent <- eta + matrix(rnorm(length(eta), 0, delta1), ncol = ncol(eta))
 # Correct zeta by adding the lambda vector
  zeta.corrected <- zeta + lambda
  return(list(Theta.lattent = Theta.lattent, zeta.corrected = zeta.corrected))
}

```


Function to create sensitivity parameter 
```{r}
create_manyLambda <- function(param, n1_catego, A, B) {
  param_data <- t(tidyr::expand_grid(param, -param))
  Lambda <- matrix(0, nrow = n1_catego - 1, ncol = ncol(param_data)) # Initialisation
  Lambda[A, ] <- param_data[1, ]
  Lambda[B-1, ] <- param_data[2, ]
  Lambda <- as.data.frame(Lambda)
  colnames(Lambda) <- paste0("lambda", 1:ncol(Lambda))
  return(Lambda)
}
manyLambda <- create_manyLambda(param = seq(0, 4, 1),n1_catego = n1_catego, A = A, B = B)

manyLambda <- manyLambda[,c(1:2,6,7:8,13)]
manyLambda #
```



```{r}
cl <- makeCluster(round(detectCores()/2))

clusterExport(cl, c("manyLambda", "BUILD","out.regression","nsim")) 
out.BUILD <- vector("list", nsim)
out.BUILD <- parLapply(cl, 1:nsim, function(s) {
  lapply(1:ncol(manyLambda), function(l) {
    BUILD(eta = out.regression[[s]]$eta, 
          zeta = out.regression[[s]]$zeta,
          lambda = manyLambda[, l], 
          seed = 10000)
  })
})
stopCluster(cl)
```





## `Étape 6` de l'algorithme


```{r}
USE <- function(data.imp,THETA, ZETA ){
  mnar <- matrix(NA, nrow = nrow(THETA), ncol = ncol(THETA)) ## Initialisation 
  data.mnar <- data.imp
  for(i in 1:ncol(THETA)){
    # Fill missing values in mnar based on THETA and the intervals defined by ZETA.
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
     # Replace missing values in x1.mnar with the corresponding values in mnar.
    data.mnar[[i]]$x1.mnar[is.na(data.mnar[[i]]$x1.mis)] <- mnar[, i]
    data.mnar[[i]]$x1.mnar <- ordered(data.mnar[[i]]$x1.mnar)
  }
  return(data.mnar)
}
```



```{r}
parallel.use <- function(s,M, manyLambda) {
  
 out.USE.chunk <- PROPORTION.chunk <- vector("list", ncol(manyLambda)) #  Initialisation 
  for (l in 1: ncol(manyLambda)) {
  # Call the USE function to generate MNAR data with the values of THETA and ZETA.
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
cl <- makeCluster(round(detectCores()/2)) 
clusterExport(cl, c("datas", "USE", "out.BUILD", "parallel.use", "out.regression", "manyLambda", "M"))
use_proportion <- parLapply(cl, 1:nsim, function(s){
   parallel.use(s,M = M, manyLambda = manyLambda)
})
stopCluster(cl)
```

```{r}
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
  mnar.table <- mar.table <- full.table <- vector("list", ncol(manyLambda)) 
  for (l in 1: ncol(manyLambda)) {
    mnar.table[[l]] <- mar.table[[l]] <- matrix(NA, nrow = M, ncol = length(unique(datas[[1]]$x1)), dimnames = list(paste0("m", 1:M), paste0("X1", 1:length(unique(datas[[1]]$x1)))))
  }
  # Compute the proportions of MNAR and MAR for each set of imputations m.
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
cl <- makeCluster(round(detectCores()/2)) 
clusterExport(cl, c("nsim", "datas", "manyLambda", "M", "PROPORTION", "parallel_proportion"))
suppressWarnings(mar.mnar.result <- parLapply(cl, 1:nsim, function(s){
  parallel_proportion(s , manyLambda = manyLambda)
}))
stopCluster(cl)
```


```{r}
mar.prop.mean <- data.frame("X1" = mar.mnar.result[[1]]$mar[[1]][, 1], Reduce("+", lapply(mar.mnar.result, function(result) result$mar[[1]][, -1])) / nsim)
mar.prop.mean <- cbind(mar.prop.mean[, 1:2], "all.m" = rowMeans(mar.prop.mean[, -c(1:2)]))
mnar.prop.mean <- lapply(1:ncol(manyLambda), function(l) {
  mnar_df <- data.frame("X1" = mar.mnar.result[[1]]$mnar[[l]][, 1], Reduce("+", lapply(mar.mnar.result, function(result) result$mnar[[l]][, -1])) / nsim)
  mnar_df <- cbind(mnar_df[, 1:2], "all.m" = rowMeans(mnar_df[, -c(1:2)]))
  return(mnar_df)
})
```




The purpose of this function is to streamline the process of loading multiple R packages by providing a convenient way to load them all at once.
```{r}
PACKAGE <- function(package){
  for(p in package){
    library(p, character.only = TRUE)
  }
}

```



```{r}
couleur <- c("red", "black",  "orange","green", "purple", "pink", "cyan","magenta", "dimgrey", "blue", "darkred",  "darkgreen", "darkblue", "slategray3","springgreen4","navy", "gold","#FF4D00", "#FF9900","#CCFF00","#3300FF" ,"#00B3FF","turquoise4", "tan","slategray1","sienna3","rosybrown4","peachpuff","mistyrose4","lightyellow","firebrick2")

cl <- makeCluster(round(detectCores()/2))
invisible(clusterCall(cl, PACKAGE, c("ggplot2", "ggthemes")))
clusterExport(cl, c("mar.prop.mean", "mnar.prop.mean", "manyLambda", "M", "nsim","couleur"))
comparaison <- parLapply(cl, 1:ncol(manyLambda), function(l) {
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

display multiple plots in a grid in R

```{r}
grid.arrange(grobs = comparaison, ncol = 2)
```




# Models d'analyse finale `regression logistique`




## Full data `Échantillons sans valeurs manquantes`

```{r,warning=FALSE, message=FALSE}
cl <- makeCluster(round(detectCores()/2)) #
clusterExport(cl, c("datas"))
#  GLM
Full.mod <- parLapply(cl, datas, function(data) {
  glm(y ~ x1 + x2 + x3 + x4, family = binomial(link = "logit"), data = data)
})
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

 Function to calculate the squared estimation variance and Function to check if a given value is within a confidence interval
```{r}
var_estime <- function(estime_ech, estime_mean) {
  # Calculate and return the squared difference between the estimated value and the mean estimate
  return((estime_ech - estime_mean)^2)
}

coverage <- function(value, CI.low, CI.upper) {
  # Use ifelse to check if the value falls within the confidence interval
  # If it does, return 1; otherwise, return 0
  return(ifelse(CI.low <= value & CI.upper >= value, 1, 0))
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

Estimate of full data

```{r}
estisimule <- full_mcar_mar(true = betas, nsim =nsim,  estimate = Full.est)
estisimule <- round(estisimule,2)
estisimule
```


##  Completes cases analysis `Mécanisme MCAR`


```{r}
cl <- makeCluster(round(detectCores()/2)) # nombre de coeurs alloués
clusterExport(cl, c("datas"))
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
stopCluster(cl)
```


Estaimate of MCAR

```{r}
estimcar <- full_mcar_mar(true = betas,nsim =nsim, estimate = MCAR.est)
estimcar <- round(estimcar,2)
estimcar
```

## Le module mice `Mécanisme MAR`

```{r}
calcul_mar <- function(data) {
    model <- with(data, glm(y ~ x1.imp + x2 + x3 + x4, family = binomial(link = "logit")))
    model_summary <- summary(mice::pool(model), conf.int = TRUE)
    est <- data.frame(
      cbind(
        model_summary[, c("estimate", "2.5 %", "97.5 %")],
        "Largeur.IC" = model_summary[, "97.5 %"] - model_summary[, "2.5 %"]
      )
    )
    
    colnames(est) <- c("estime",  "Borne.inf", "Borne.sup","Largeur.IC")
    rownames(est) <-  c("Intercept", paste0("X1", 2:6), "X2", "X32", "X33", "X41")
    MAR.est <- est
    return(MAR.est = MAR.est)
  }
```


```{r}
cl <- makeCluster(round(detectCores()/2)) # nombre de coeur alloué
clusterExport(cl, c("mice.imp","calcul_mar"))
invisible(clusterCall(cl, PACKAGE, c("mice")))
tout_MAR.est <- parLapply(cl,  1:nsim, function(s){
  calcul_mar(data = mice.imp[[s]])
})
stopCluster(cl)
```
Estimate of MAR
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
cl <- makeCluster(round(detectCores()/2))
clusterExport(cl, c("nsim", "manyLambda", "out.USE",  "prara.mnar"))
invisible(clusterCall(cl, PACKAGE, c("mice", "mitml")))

tout_MNAR.est <- parLapply(cl, 1:nsim, function(s){
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


Estimate of MNAR

```{r}
estimnar <- lapply(1: ncol(manyLambda), function(l) {
  round(mnar(true = betas, estimate = tout_MNAR.est, l = l, nsim = nsim),2)
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
graph1 <- estime_graph(data = estisimule, title = " Beta True  vs Beta Simule ", labels = c( "True", "Simule"),colors = couleur[1:2])
graph2 <- estime_graph(data = estimcar, title = " Beta True  vs Beta Mcar ", labels = c( "True", "Mcar"),colors = couleur[c(1,3)])
graph3 <- estime_graph(data = estimar, title = " Beta True  vs Beta Mar ", labels = c( "True", "Mar"),colors = couleur[c(1,4)])
graph <- list(graph1,graph2,graph3)


graph.mnar <- lapply(1:ncol(manyLambda), function(l) {
  estime_graph(data = estimnar[[l]], title = paste("Beta True vs",paste( "Beta", paste0("Mnar",l))), labels = c("True", paste0("Mnar",l)), colors = c(couleur[1], couleur[5:(ncol(manyLambda)+4)][l]))
})

```

```{r}
grid.arrange(grobs = c(graph, graph.mnar), ncol = 3)
```


```{r}
all_estime_biais_sd_coverage <- function(betas, estisimule, estimcar, estimar, estimnar) {
  
  # Resumé des estimés 
  all_estimate <- data.frame(
     betas,
    estisimule["estime"],
      estimcar["estime"],
       estimar["estime"],
       sapply(1:length(estimnar),function(l) estimnar[[l]]["estime"])
  )
  
  # Variabilite
  all_sd_estime <- data.frame(
    estisimule["sd_estime"],
      estimcar["sd_estime"],
       estimar["sd_estime"],
       sapply(1:length(estimnar),function(l) estimnar[[l]]["sd_estime"])
  )

  #  Résumé du biais relatif
  all_relative_biais <- as.data.frame(matrix(NA, nrow = nrow( all_estimate), ncol = ncol( all_estimate)))
  for (i in 1:ncol( all_estimate)) {
   all_relative_biais[[i]] <- round(100 * (( all_estimate[, i] -  all_estimate[, 1]) /  all_estimate[, 1]),2)
  }
  
  
   # Proportion des ICs qui contiennent la valeur
  all_coverage <- data.frame(
    estisimule["couvert"],
      estimcar["couvert"],
       estimar["couvert"],
       sapply(1:length(estimnar),function(l) estimnar[[l]]["couvert"])
  )
  
   colnames(all_estimate) <- colnames(all_relative_biais) <-  c("True", "Simule","Mcar","Mar", paste0("Mnar", 1:length(estimnar)))
   colnames(all_sd_estime) <- colnames(all_coverage) <- colnames(all_estimate)[-1]
  rownames(all_estimate) <- rownames( all_sd_estime)<- row.names(all_relative_biais) <- rownames(estisimule)
  
  return(list(all_estimate = all_estimate, all_sd_estime = all_sd_estime, all_relative_biais = all_relative_biais, all_coverage = all_coverage))
}



all_result <- all_estime_biais_sd_coverage(betas = betas, estisimule = estisimule, estimcar = estimcar, estimar = estimar, estimnar = estimnar)
```



all estimate
```{r}
all_estimate <- all_result$all_estimate
all_estimate
```

all relative biais
```{r}
all_relative_biais <- all_result$all_relative_biais
all_relative_biais
```


```{r}
df_biais <- data.frame("coefficiant" =rownames(all_relative_biais), all_relative_biais)
rownames(df_biais)  <- NULL

# Utiliser la fonction gather pour regrouper les données

df_biais_long <- tidyr::gather(df_biais, Method, Bias, -coefficiant)

labels_methods <- unique(df_biais_long$Method)

# Créer le graphique ggplot avec des lignes connectant les points de même couleur
gg <- ggplot(df_biais_long, aes(x = coefficiant, y = Bias, color = Method, group = Method)) +
  geom_point() +
  geom_line() +
  scale_color_manual(
    values = couleur[1:length(labels_methods)],
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


# Afficher le graphique
print(gg)
```

all SE 

```{r}
## Sd  des estimés
all_sd_estime <- all_result$all_sd_estime
all_sd_estime
```




## all coverage

```{r}
## Coverage
all_coverage <- all_result$all_coverage
all_coverage
```



