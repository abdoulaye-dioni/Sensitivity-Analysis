Imputation is a widely recognized method for addressing missing data. Among the packages used for this in hierarchical contexts, jomo is important. However, when working with real data, comparing missing data under the MCAR (Missing Completely at Random) and MAR (Missing at Random) assumptions might not accurately capture the true nature of the missing data mechanism. Therefore, it is useful to consider the MNAR (Missing Not at Random) case using a sensitivity analysis approach to evaluate deviation from MAR because, in some situations, the jomo package can produce potentially biased estimates.
The objective is to develop a sensitivity analysis method for ordinal variables where certain overrepresented categories are potentially missing, possibly depending on the response variable. Additionally, the proportions of missing values in each category of an ordinal variable are likely to vary across different strata of the categorical exposure variable of interest.
The proposed method is based on an algorithm presented in the article titled "Sensitivity analysis method in the presence of a missing not at random ordinal independent variable."
The method should be used with caution, as it heavily relies on the assumption that certain specific (overrepresented) categories of an ordinal variable tend to be missing in each stratum of the categorical exposure variable of interest.


```{r}
#sink("hierachical.txt")
ti <- Sys.time()
```


```{r}
library(parallel)
library(mitml)
library(lme4)
library(ggplot2)
library(gridExtra)
library(ggthemes)
library(jomo)
```





# Simulation

The `simulation_hierarchical` function generates a set of simulated data with hierarchical structures from a generalized linear mixed-effects model. Without loss of generality, we considered the case of a random intercept.

```{r}
simulation_hierachical <- function(length_x1, length_x2, formula_x1, formula_y, para_x1,  betas,n_clus,n_obs,sd_U) {
  x2 <- factor(sample(1:length_x2, n_clus*n_obs, replace = TRUE, prob = rep(0.25,4)))
  eta_x1 <- model.matrix(as.formula(formula_x1), data.frame( x2 = x2))[,-1] %*% as.matrix(para_x1)
  prob_x1 <-  t(apply(eta_x1, 1, function(x) exp(x - max(x)) / sum(exp(x - max(x)))))
  x1 <- ordered(apply(prob_x1, 1, function(p) sample(1:length(p), 1, prob = p)))
  clus <- rep(1:n_clus, each=n_obs)
  U <- rnorm(n_clus, mean=0, sd = sd_U)
  ZU <- U[clus]
  eta_y <- model.matrix(as.formula(formula_y), data = data.frame(x1 = x1, x2 = x2)) %*%betas + ZU
  prob_y <- plogis(eta_y)
  y <- as.factor(rbinom(n_clus*n_obs, 1, prob_y))
  
  return(data.frame(id = seq_len(n_clus*n_obs), y = y, x1 = x1,  x2 = x2, clus = clus))
}
```







## Parameters 

* `length_x1`: the number of levels for the ordinal variable x1.
* `length_x2`: the number of levels for the categorical variable x2.
* `formula_x1`: the formula used to generate the values for x1.
* `formula_y`: the formula used to generate the values for y.
* `para_x1`: a matrix of parameters used to generate x1.
* `betas`: a vector of coefficients used to generate y.
* `n_clus`: the number of clusters or groups.
* `n_obs`: the number of observations per cluster.
* `sd_U`: the standard deviation for the random effects of the clusters.

```{r}
set.seed(100)
nsim = 500
n_clus = 10
n_obs = 200
length_x1 <- 3
length_x2 = 4
sd_U = 0.45
A = 1
B = 3


formula_x1 <- "~ x2"
para_x1 <- matrix(c(rep(2,3), rep(1,3),rep(2.5,3) ),ncol = 3)


formula_y <- "~ x1 + x2"
betas <- c(-1, 1, -2, 2,  1, 2)
```



```{r}
datas <- lapply(1:nsim, function(simulation) {
    data <- simulation_hierachical(length_x1, length_x2, formula_x1, formula_y, para_x1,  betas,n_clus,n_obs,sd_U)
  data$x1.mis <- data$x1
  return(data)
  })
```




##  Missing Data on the extreme categories of ordinal Variable (MNAR)

Creation of probabilities for missing values to be applied to the extreme categories of the ordinal variable for each stratum of the categorical exposure variable of interest.

```{r}
missing_prob <- data.frame(matrix(c(0.2,0.3,0.1,0.4,0.4,0.1,0.1,0.3), nrow = 2, byrow = FALSE,dimnames = list(c("A","B"),paste0("proba",1:length_x2))))
missing_prob 
```

The `simulate_mnar` function simulates Missing Not At Random (MNAR) values on the ordinal variable for each stratum of the categorical exposure variable of interest.


```{r}
simulate_mnar <-  function(data, proba, cat_var) {
  
  mnar_by_strata <- function(data, repY, id, ord_var,probA, probB, A, B, seed) {
  set.seed(seed)
  id.A <- data[[id]][data[[repY]] %in% 1 & data[[ord_var]] %in% A]
  sample.A <- sample(id.A, size = round(length(id.A) * probA), replace = FALSE)
  data[[ord_var]][data[[id]]%in%sample.A] <- NA
  id.B <- data[[id]][data[[repY]] %in% 0 & data[[ord_var]] %in% B]
  sample.B <- sample(id.B, size = round(length(id.B) * probB), replace = FALSE)
  data[[ord_var]][data[[id]]%in%sample.B] <- NA
  return(data)
}
  result <- lapply(unique(data[[cat_var]]), function(i) {
    mnar_by_strata(data = subset(data, subset = data[[cat_var]] == i), 
                  repY = "y", id = "id", ord_var = "x1.mis", 
                  probA = proba[[i]][1], probB = proba[[i]][2], 
                  A = A, B = B, seed = 100)
  })
  return(do.call(rbind, result))
}
```




```{r}
for(s in 1:nsim){
  resultat_s <- simulate_mnar(data = datas[[s]], cat_var = "x2", proba = missing_prob)
  datas[[s]] <- resultat_s[order(resultat_s$id),]
  datas[[s]]$x1.imp <- datas[[s]]$x1.mis
}
```




# Multiple Imputation with the jomo package


This `PACKAGE` function loads multiple R packages using their names provided in a list or vector for parallel computation in R.

```{r}
PACKAGE <- function(package){
  for(p in package){
    library(p, character.only = TRUE)
  }
}
```





```{r}
M <- 15
nburn <- 1000
nbetween <- 1000
cl <- makeCluster(round(detectCores()/2))
formula = y ~ x1.imp + x2 + (1|clus) 
clusterExport(cl, c("datas","M", "nburn", "nbetween", "formula","nsim"))
invisible(clusterCall(cl, PACKAGE, c("jomo","mitml")))

# Calcule en parallele
mar_data <- parLapply(cl, 1:nsim, function(s) {
  jomo.glmer(
    formula = as.formula(formula),
    data = datas[[s]][, c("y", "x1.imp", "x2", "clus")],
    nburn = nburn,
    nbetween = nbetween,
    nimp = M,
    output = 0
  )
})

mar_data <- parLapply(cl, mar_data, function(jomo_imp) {
  jomo2mitml.list(jomo_imp)
})

stopCluster(cl) 
```


## Sensitivity Analysis

Creating a database to contain the data modified under MNAR, which will be updated throughout the code.

```{r}
mnar_data <- lapply(1:nsim, function(s) {
  lapply(1:M, function(m) {
    cbind(
      mar_data[[s]][[m]],
      datas[[s]][, c("x1", "x1.mis")]
    )
  })
})

```

## Ordinal Regression

The `ordinal_regression` function performs an ordinal regression with a probit link function on a given dataset and generates results that include the intercepts ("zeta") and linear predictors ("eta"). The function adds Gaussian noise to the linear predictors ("eta"). Indeed, the coefficients `beta` and the intercepts `zeta` are potentially biased and imprecise. The addition of Gaussian noise accounts for the imprecision of data imputed under MAR.

```{r}
ordinal_regression <- function(data,M, length_x1, seed = NULL) {
  set.seed(seed)
  ord.regression <- lapply(data, function(imp_data) {
    MASS::polr(x1.imp ~ y + x2, method = "probit", Hess = TRUE, data = imp_data)
  })
  
  zeta <- sapply(ord.regression, function(ord_reg) {
    as.numeric(ord_reg$zeta)
  })
  
  data <- lapply(seq_along(data), function(m) {
    data[[m]]$eta <- as.numeric(ord.regression[[m]]$lp + rnorm(nrow(data[[1]]),0,1.5))
    return(data[[m]])
  })
  dimnames(zeta) <- list(paste0("k", 1:(length_x1 -1)), paste0("m", 1:M))
  
  return(list(data = data, zeta = zeta))
}
```




```{r}
cl <- makeCluster(round(detectCores()/2))
clusterExport(cl, c("ordinal_regression", "mnar_data", "M", "length_x1")) 
out_ordinal_regression <- parLapply(cl, 1:nsim, function(s) {
  ordinal_regression(
   data = mnar_data[[s]],
   M=M,
   seed = 123,
   length_x1 =length_x1
  )
})

stopCluster(cl) 
```








## Modification of the Zeta Intercepts

Choice of sensitivity parameter vectors that will be used to adjust the linear predictor, taking into account additional information that the categories of the ordinal variable tend to be missing.


```{r}
manydelta <- data.frame( delta1 = c(0.5,0), delta2 = c(0,-0.5), delta3 = c(0,-1.5),delta4 = c(0,-2))
manydelta 
```



```{r}
compute_zeta_new <- function(zeta_old, manydelta) {
  zeta_new <- lapply(seq_along(manydelta), function(l) {
    zeta_old + manydelta[, l]
  })
  return(zeta_new)
}
```



```{r}
cl <- makeCluster(round(detectCores() /2))
clusterExport(cl, c("compute_zeta_new", "out_ordinal_regression", "manydelta", "M"))
zeta.new <- parLapply(cl, 1:nsim, function(s) {
    compute_zeta_new(out_ordinal_regression[[s]]$zeta, manydelta)
})

stopCluster(cl)
```




## Modification of Data Imputed under MAR

Updating the `mnar_data` database.



```{r}
mnar_data <- lapply(1:nsim, function(s) {
  lapply(1:M, function(m) {
    mnar_columns <- out_ordinal_regression[[s]]$data[[m]][, c(rep("x1.mis", ncol(manydelta) +3)), drop = FALSE]
    colnames(mnar_columns) <- c(paste0("delta", 1:ncol(manydelta)),paste0("mnar", 1:3))
    regression_data <- cbind(out_ordinal_regression[[s]]$data[[m]], mnar_columns)
    return(regression_data)
  })
})
```






The `compute_strata` function creates strata from the categorical exposure variable of interest.



```{r}
compute_strata <- function(data, var) {
  strata <- list()
  for (j in levels(data[[var]])) {
    strata[[j]] <- subset(data, data[[var]] == j) 
  }
  return(strata)  
}
```


```{r}
cl <- makeCluster(round(detectCores()/2))
clusterExport(cl, c( "mnar_data", "compute_strata"))
my_strata <- parLapply(cl, 1:nsim, function(s) {
      lapply(mnar_data[[s]], function(dat) {
  compute_strata(data = dat, var = "x2")})
})

stopCluster(cl)

```






The `modified_mar` function modifies a dataset imputed under MAR by considering the new `eta and zeta` values obtained earlier. This creates a new dataset imputed under MNAR and is consistent with the information that certain specific categories of the ordinal variable tend to be missing. The function is applied to each stratum of the categorical exposure variable to account for the variability in the proportions of missing data within each stratum

```{r}
modified_mar <- function(data, eta, zeta, m,seed = NULL) {
  set.seed(seed)
  for (l in 1:ncol(manydelta)) {
    for (i in 1:length(eta)) {
      if(is.na(data[i, paste0("delta", l)])) {
data[, paste0("delta", l)][i][eta[i] <= zeta[[l]][,m][1]] <- 1
data[, paste0("delta", l)][i][eta[i] > zeta[[l]][,m][1] & eta[i] <= zeta[[l]][,m][2]] <- 2
data[, paste0("delta", l)][i][eta[i] > zeta[[l]][,m][2]] <- 3
    }
  }
  data[, paste0("delta", l)][is.na(data$x1.mis)] <- sample(data[, paste0("delta", l)][is.na(data$x1.mis)])
  }
  return(data)
}
```



```{r}
cl <- makeCluster(round(detectCores() /2))
clusterExport(cl, c("my_strata", "modified_mar", "zeta.new","manydelta", "M", "length_x2"))
mnar_data <- parLapply(cl, 1:nsim, function(s) {
  lapply(1:M, function(m) {
    lapply(1:length_x2, function(a) {
      modified_mar(data = my_strata[[s]][[m]][[a]],
                    eta = my_strata[[s]][[m]][[a]]$eta,
                    zeta = zeta.new[[s]],
                    m=m,
                   seed = 123)
    })
  })
})

stopCluster(cl)
```



## Diagnosing Data Imputed under MAR vs. MNAR for each stratum.

The `calculate_missing_proportions` function calculates and returns the proportions of categories imputed under `mar`, modified under `mnar`, and in the case where there are no missing categories (`full data`) for each stratum of the categorical exposure variable of interest.

```{r}
calculate_missing_proportions <- function(data, nsim, manydelta,length_x1, length_x2) {
  #  MAR
  proportion_mar <- function(s, M) {
    matrices <- lapply(1:M, function(m) {
    sapply(1:length_x2, function(a) {
     prop.table(table(data[[s]][[m]][[a]]$x1.imp[is.na(data[[s]][[m]][[a]]$x1.mis)])) * 100
    })
    })
    return(t(Reduce("+", matrices) / M))
  }
  
  #  MNAR
  
  proportion_mnar <- function(s,l,M) {
    matrices <- lapply(1:M, function(m) {
      sapply(1:length_x2, function(a) {
        prop.table(table(data[[s]][[m]][[a]][,paste0("delta", l)][is.na(data[[s]][[m]][[a]]$x1.mis)])) * 100
      })
    })
    return(t(Reduce("+", matrices) / M))
  }
  
  #  full data
  proportion_full <- function(s) {
    matrices <- sapply(1:length_x2, function(a) {
      prop.table(table(data[[s]][[1]][[a]]$x1[is.na(data[[s]][[1]][[a]]$x1.mis)])) * 100
    })
    return(t(matrices))
  }
  
  # MAR proportions
  mar <- round(Reduce("+", lapply(1:nsim, function(s) { proportion_mar(s, M) })) / nsim, 2)
  
  # MNAR proportion
  
mnar <- lapply(1:ncol(manydelta), function(l) {
    round(Reduce("+", lapply(1:nsim, function(s) { proportion_mnar(s, l, M) })) / nsim, 2)
  })

  
  # full data proportions 
  full <- round((Reduce("+", lapply(1:nsim, function(s) { proportion_full(s) })) / nsim), 2)
  
  #  final table
  my.table <- lapply(1:length_x2, function(a) {
    matrix(NA, nrow = length_x1, ncol = 2+ ncol(manydelta),
           dimnames = list(paste0("X1", 1:length_x1), c("full", "mar", paste0("delta", 1:ncol(manydelta)))))
  })
  
  for (a in 1:length_x2) {
    my.table[[a]][,1] <- full[a,]
    my.table[[a]][,2] <- mar[a, ]
    
    for (l in 1:ncol(manydelta)) {
      my.table[[a]][,2 + l] <- mnar[[l]][a,]
    }
  }
  
  return(my.table)
}
```




```{r}
my.table <- calculate_missing_proportions(mnar_data, nsim, manydelta,length_x1, length_x2)

my.table
```



The `generate_plot` function creates a graphical representation of these different proportions for each stratum of the categorical exposure variable of interest.

```{r}
couleur <- c("red", "black",  "blue",  "orange","green", "magenta","darkgreen", "purple", "cyan", "#CCFF00")
# Fonction pour les graphes
generate_plot <- function(data, l, category_title, delta_value) {
  ggplot(data) +
    geom_path(aes(x = X1, y = full, group = 1, col = "black")) +
    geom_point(aes(x = X1, y = full, col = "black")) +
    geom_path(aes(x = X1, y = mar, group = 1, col = "orange")) +
    geom_point(aes(x = X1, y = mar, col = "orange")) +
    geom_path(aes(x = X1, y = !!sym(paste0("delta", l)), group = 1, col = "magenta")) +
    geom_point(aes(x = X1, y = !!sym(paste0("delta", l))), col ="magenta") +
    labs(
      y = "Proportion",
      color = "Legend"
    ) +
    scale_color_manual(
      name = "Legend", 
      values = c("black", "orange", "magenta"), 
      breaks = c("black", "orange", "magenta"), 
      labels = c("Simule", "Mar", paste0("Delta", l))
    ) +
    theme_clean(base_size = 6) +
    ggtitle(bquote(.(category_title) * "  ;  " * delta[.(l)] == .(delta_value))) + 
    theme(plot.title = element_text(hjust = 0.5)) +
     theme(axis.text.x = element_text(angle = 45, hjust = 0.5)) +
    ylim(0, 100)
}

```





```{r}
resulta <- list()

for (a in 1:length(my.table)) {
  strata_graph <- data.frame(X1 = row.names(my.table[[a]]), my.table[[a]])
  row.names(strata_graph) <- NULL
  resulta[[a]] <- lapply(1:ncol(manydelta), function(l) generate_plot(data = strata_graph, l = l, category_title = paste("x2=", a), delta_value = paste(as.vector(manydelta[, l]), collapse = ", ")))
}
```


```{r}
grid.arrange(grobs = do.call(c, lapply(1:length_x2, function(i) eval(parse(text = paste0("resulta[[", i, "]]"))))), ncol = ncol(manydelta))
```


Different combinations of data modified under MNAR to form `mna1 and mnar2` that are plausible for each stratum of the categorical exposure variable of interest.

```{r}
for(s in 1:nsim){
  for(m in 1:M){
    mnar_data[[s]][[m]][[1]][["mnar1"]] <- mnar_data[[s]][[m]][[1]][["delta4"]]
    mnar_data[[s]][[m]][[2]][["mnar1"]] <- mnar_data[[s]][[m]][[2]][["delta3"]]
    mnar_data[[s]][[m]][[3]][["mnar1"]] <- mnar_data[[s]][[m]][[3]][["delta1"]]
    mnar_data[[s]][[m]][[4]][["mnar1"]] <- mnar_data[[s]][[m]][[4]][["delta3"]]
  }
}


for(s in 1:nsim){
  for(m in 1:M){
    mnar_data[[s]][[m]][[1]][["mnar2"]] <- mnar_data[[s]][[m]][[1]][["delta2"]]
    mnar_data[[s]][[m]][[2]][["mnar2"]] <- mnar_data[[s]][[m]][[2]][["delta4"]]
    mnar_data[[s]][[m]][[3]][["mnar2"]] <- mnar_data[[s]][[m]][[3]][["delta1"]]
    mnar_data[[s]][[m]][[4]][["mnar2"]] <- mnar_data[[s]][[m]][[4]][["delta4"]]
  }
}

for(s in 1:nsim){
  for(m in 1:M){
    mnar_data[[s]][[m]][[1]][["mnar3"]] <- mnar_data[[s]][[m]][[1]][["delta4"]]
    mnar_data[[s]][[m]][[2]][["mnar3"]] <- mnar_data[[s]][[m]][[2]][["delta4"]]
    mnar_data[[s]][[m]][[3]][["mnar3"]] <- mnar_data[[s]][[m]][[3]][["delta4"]]
    mnar_data[[s]][[m]][[4]][["mnar3"]] <- mnar_data[[s]][[m]][[4]][["delta4"]]
  }
}
```


```{r}
mnar_data <- lapply(lapply(1:nsim, function(s) {
  lapply(1:M, function(m) {
    do.call(rbind, mnar_data[[s]][[m]][1:length_x2])
  })
}), function(s) {
  lapply(s, function(my_ord) {
    my_ord[order(my_ord$id), ]
  })
})
```




# Final Analysis 

## Full Data 

generalized linear mixed model on simulated data without missing values.



```{r}
control.param =glmerControl(optimizer = "Nelder_Mead", 
optCtrl = list(maxfun = 100000), check.conv.grad = .makeCC("warning",
tol = 7e-2, relTol = NULL), check.conv.hess =.makeCC(action = "warning", tol = 1e-10),
boundary.tol = 1e-10)
```


```{r}
cl <- makeCluster(round(detectCores()/2)) 
clusterExport(cl, c("datas", "control.param", "length_x1", "length_x2"))
invisible(clusterCall(cl,PACKAGE, "lme4"))
#  GLMM 
Full.mod <- parLapply(cl, datas, function(data) {
  glmer(y ~ x1 + x2 + (1 | clus), family = binomial, data = data, control = control.param)
})
# Estimate and CI
Full.est <- parLapply(cl, Full.mod, function(model) {
  estime <- coef(summary(model))[,"Estimate"]
  conf_intervals <- confint.merMod(model, method = "Wald")[-1,]
  est_df <- data.frame(
    "estime" = estime,
    "Borne.inf" = conf_intervals[, "2.5 %"],
    "Borne.sup" = conf_intervals[, "97.5 %"],
    "Largeur.IC" = conf_intervals[, "97.5 %"] - conf_intervals[, "2.5 %"]
  )
  rownames(est_df) <- c("Intcept", paste0("X1", 2:length_x1),  paste0("X2", 2:length_x2))
  return(est_df)
})
stopCluster(cl)
```


The `var_estime` function calculates the empirical standard deviation, and the `coverage` function calculates the coverage proportions of confidence intervals.

```{r}
var_estime <- function(estime_ech, estime_mean){
  return((estime_ech - estime_mean)^2)
}

coverage <- function(value, CI.low, CI.upper){
  ifelse(CI.low <= value & CI.upper >= value,1,0)
}
```


The `full_mcar_mar` function performs statistical calculations to evaluate the estimates, their variance, their confidence intervals, and their coverage proportions. It is applicable in the case of simulated "full data," data imputed under "MAR," and modified "MNAR."

```{r}
full_mcar_mar <- function(true, estimate, nsim){
  estime <- var_by_sim <-  couvert <- Largeur.IC <- data.frame(matrix(NA, nrow = nrow(estimate[[1]]),ncol = nsim))
  for(s in 1:nsim){
    estime[,s] <- estimate[[s]]["estime"]
    Largeur.IC[,s] <- estimate[[s]]["Largeur.IC"]
  var_by_sim[,s] <- var_estime(estime_ech = estimate[[s]][,"estime"], estime_mean = rowMeans(sapply(estimate, function(x) x[,"estime"])) )
  couvert[,s] <- coverage(value = true, CI.low = estimate[[s]]["Borne.inf"], CI.upper = estimate[[s]]["Borne.sup"])
  }

  estime_IC <- data.frame(
    "true" = true,
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


```{r}
estisimule <- full_mcar_mar(true = betas, nsim =nsim,  estimate = Full.est)
round(estisimule,2)
```

## MCAR Mechanism

generalized linear mixed model under MCAR

```{r}
cl <- makeCluster(round(detectCores()/2)) 
clusterExport(cl, c("datas", "control.param", "length_x1", "length_x2"))
invisible(clusterCall(cl,PACKAGE, "lme4"))
#  GLMM
MCAR.mod <- parLapply(cl, datas, function(data) {
  glmer(y ~ x1.mis + x2 + (1 | clus), family = binomial, na.action = na.omit, data = data, control = control.param)
})

MCAR.est <- parLapply(cl, MCAR.mod, function(model) {
  estime <- coef(summary(model))[,"Estimate"]
  conf_intervals <- confint.merMod(model, method = "Wald")[-1,]
  est_df <- data.frame(
    "estime" = estime,
    "Borne.inf" = conf_intervals[, "2.5 %"],
    "Borne.sup" = conf_intervals[, "97.5 %"],
    "Largeur.IC" = conf_intervals[, "97.5 %"] - conf_intervals[, "2.5 %"]
  )
  rownames(est_df) <- c("Intcept", paste0("X1", 2:length_x1),  paste0("X2", 2:length_x2))
  return(est_df)
})

stopCluster(cl)

```




```{r}
estimcar <- full_mcar_mar(true = betas,nsim =nsim, estimate = MCAR.est)
round(estimcar,2)
```



## MAR Mechanism

generalized linear mixed model under MAR

```{r}
calcul_mar <- function(data,length_x1,length_x2) {
    MAR.mod <- with(data, glmer(y ~ x1.imp + x2 + (1 | clus), family = binomial, control = control.param))
    estime <- testEstimates(MAR.mod, extra.pars = TRUE)$estimates[,"Estimate"]
    conf_intervals <- confint(testEstimates(MAR.mod, extra.pars = TRUE))
    MAR.est <- data.frame(
      "estime" = estime,
    "Borne.inf" = conf_intervals[, "2.5 %"],
    "Borne.sup" = conf_intervals[, "97.5 %"],
    "Largeur.IC" = conf_intervals[, "97.5 %"] - conf_intervals[, "2.5 %"]
    )
    rownames(MAR.est) <- c("Intcept", paste0("X1", 2:length_x1),  paste0("X2", 2:length_x2))
  return(MAR.est = MAR.est)
}
```


```{r}
cl <- makeCluster(round(detectCores()/2))
invisible(clusterCall(cl, PACKAGE, c("mitml", "lme4")))
clusterExport(cl, c("calcul_mar","mar_data", "control.param", "length_x1", "length_x2"))

tout_MAR.est <- parLapply(cl, 1:nsim, function(s){
  calcul_mar(mar_data[[s]],length_x1,length_x2)})
stopCluster(cl)
```



```{r}
estimar <- full_mcar_mar(true = betas,nsim = nsim, estimate = tout_MAR.est)
round(estimar,2)
```




## MNAR Mechanism

generalized linear mixed model under MNAR

```{r}
mnar_data <-   lapply(1:nsim, function(s){mitml::as.mitml.list(mnar_data[[s]])})
```


```{r}
calcul.mnar <- function(data, manydelta, length_x1, length_x2) {
  #variables <- c(paste0("delta", 1:ncol(manydelta)), paste0("mnar", 1:2))
   variables <- paste0("mnar", 1:3)
  MNAR.mod <- lapply(variables, function(variable) {
    with(data, glmer(as.formula(paste("y ~", variable, "+ x2 + (1 | clus)")), family = binomial, control = control.param))
  })
  
  MNAR.est <- lapply(MNAR.mod, function(mod) {
    estime <- testEstimates(mod, extra.pars = TRUE)$estimates
    conf_intervals <- confint(testEstimates(mod, extra.pars = TRUE))
    
    data.frame(
      "estime" = estime[,"Estimate"],
      "Borne.inf" = conf_intervals[, "2.5 %"],
      "Borne.sup" = conf_intervals[, "97.5 %"],
      "Largeur.IC" = conf_intervals[, "97.5 %"] - conf_intervals[, "2.5 %"]
    )
  })
  
  MNAR.est <- lapply(MNAR.est, function(est) {
    est$term <- NULL
    rownames(est) <- c("Intcept", paste0("X1", 2:length_x1),  paste0("X2", 2:length_x2))
    return(est)
  })
  
  return(MNAR.est)
}
```


```{r}
cl <- makeCluster(round(detectCores()/2))
clusterExport(cl, c( "manydelta", "mnar_data", "calcul.mnar", "control.param", "length_x1","length_x2"))
invisible(clusterCall(cl, PACKAGE, c("lme4", "mitml")))

tout.MNAR.est <- parLapply(cl, 1:nsim, function(s){
  calcul.mnar(mnar_data[[s]], manydelta, length_x1, length_x2)
})

stopCluster(cl)
```


The `mnar` function is similar to the `full_mcar_mar` function. The difference is that it focuses on estimates under `mnar`.

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




```{r}
estimnar <- lapply(1: length(tout.MNAR.est[[1]]), function(l) {
  mnar(true = betas, estimate = tout.MNAR.est, l = l, nsim = nsim)
})

lapply(1:length(tout.MNAR.est[[1]]), function(l){
  round(estimnar[[l]],2)
})
```









```{r}

tf <- Sys.time()
code_time = tf - ti
code_time
```


# Graphs of Estimates and Key Results

## Graphs of Estimates

The `estime_graph` function allows for the graphical representation of estimates and their confidence intervals in the context of "full data," "MAR," and "MNAR."


```{r}
estime_graph <- function(data, title, labels, colors) {
  p <- ggplot() +
    #geom_line(data = as.data.frame(data), aes(x = rownames(data), y = true, group = 1, col = labels[1])) +
    geom_point(data = as.data.frame(data), aes(x = rownames(data), y = true, color = labels[1], shape = labels[1]), size = 2) +
    geom_point(data = as.data.frame(data), aes(x = rownames(data), y = estime, color = labels[2], shape = labels[2]), size = 2) +
    geom_errorbar(data = as.data.frame(data), aes(x = rownames(data), y = estime, ymin = Borne.inf, ymax = Borne.sup, color = labels[2]), width = 0.2) +
    labs(title = title, x = "Coefficients", y = "Confidence intervals") +
    scale_color_manual(name = "Legend", values = setNames(colors, labels)) +
    scale_shape_manual(name = "Legend", values = setNames(c(17, 19), labels)) +
    theme_clean(base_size = 9) +
    ylim(range(data)[1], range(data)[2])+
    theme(axis.text.x = element_text(angle = 45, hjust = 0.8)) 
  return(p)
}
```



```{r}
graph1 <- estime_graph(data = estisimule, title = " Beta True  vs Beta Simule ", labels = c( "True", "Simule"),colors = couleur[1:2])
graph2 <- estime_graph(data = estimcar, title = " Beta True  vs Beta Mcar ", labels = c( "True", "Mcar"),colors = couleur[c(1,3)])
graph3 <- estime_graph(data = estimar, title = " Beta True  vs Beta Mar ", labels = c( "True", "Mar"),colors = couleur[c(1,4)])
graph <- list(graph1,graph2,graph3)


graph.mnar <- lapply(1:length(estimnar), function(l) {
  estime_graph(data = estimnar[[l]], title = paste("Beta True vs",paste( "Beta", paste0("mnar",l)  )), labels = c("True", paste0("mnar",l)), colors = c(couleur[1], couleur[5:(ncol(manydelta)+5)][l]))
})
```


```{r}
grid.arrange(grobs = c(graph, graph.mnar), ncol = length_x1)
```


## Summary of Key Results

This `all_estime_biais_sd_coverage` function summarizes the relevant statistics.

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
  
  # colnames(all_estimate) <- colnames(all_relative_biais) <-  c("True", "Simule","Mcar","Mar", paste0("delta", 1:ncol(manydelta)), paste0("mnar", 1:2))
  colnames(all_estimate) <- colnames(all_relative_biais) <-  c("True", "Simule","Mcar","Mar", paste0("mnar", 1:3))
   colnames(all_sd_estime) <- colnames(all_coverage) <- colnames(all_estimate)[-1]
  rownames(all_estimate) <- rownames( all_sd_estime)<- row.names(all_relative_biais) <- rownames(estisimule)
  
  return(list(all_estimate = all_estimate, all_sd_estime = all_sd_estime, all_relative_biais = all_relative_biais, all_coverage = all_coverage))
}

```

```{r}
all_result <- all_estime_biais_sd_coverage(betas = betas, estisimule = estisimule, estimcar = estimcar, estimar = estimar, estimnar = estimnar)
```


### Estimates 

```{r}
all_estimate <- all_result$all_estimate
round(all_estimate,2)
```

### Relative bias

```{r}
all_relative_biais <- all_result$all_relative_biais
all_relative_biais
```


### Graph of Biases

```{r}
df_biais <- data.frame("coefficiant" =rownames(all_relative_biais), all_relative_biais)
rownames(df_biais)  <- NULL

df_biais_long <- tidyr::gather(df_biais, Method, Bias, -coefficiant)

labels_methods <- unique(df_biais_long$Method)

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
  ylim(-60,60)
```



### Graph of Biases 

```{r}
print(gg)
```

### Standard Deviation of the Estimates 

```{r}
all_sd_estime <- all_result$all_sd_estime
round(all_sd_estime,2)
```



### Coverage 

```{r}
all_coverage <- all_result$all_coverage
round(all_coverage,2)


#sink()
```
